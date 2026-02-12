#!/usr/bin/env python3
"""
CrossCell Benchmark Output Validator

验证转换输出文件的正确性：细胞数、基因数、矩阵维度。
读取工具的结果 JSON，对每个成功的转换检查输出文件。

用法:
  python3 /benchmark/scripts/validate_output.py --tool crosscell
  python3 /benchmark/scripts/validate_output.py --tool all
"""

import argparse
import json
import os
import subprocess
import sys

RESULTS_DIR = "/benchmark/results"
CONFIG_DIR = "/benchmark/config"
MATRIX_FILE = os.path.join(CONFIG_DIR, "test_matrix.json")
OUTPUT_TMP_DIR = "/tmp/benchmark_output"


def load_matrix(path):
    with open(path) as f:
        return json.load(f)


def load_results(tool_name):
    path = os.path.join(RESULTS_DIR, f"{tool_name}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def validate_h5ad(filepath, expected_cells, expected_genes):
    """用 anndata 验证 H5AD 文件"""
    try:
        script = f"""
import anndata
adata = anndata.read_h5ad("{filepath}")
print(f"{{adata.n_obs}},{{adata.n_vars}}")
"""
        result = subprocess.run(
            ["python3", "-c", script],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            return {"valid": False, "error": result.stderr[:300]}

        line = result.stdout.strip()
        n_obs, n_vars = line.split(",")
        n_obs, n_vars = int(n_obs), int(n_vars)

        ok = True
        issues = []
        if expected_cells and n_obs != expected_cells:
            issues.append(f"cells: expected {expected_cells}, got {n_obs}")
            ok = False
        if expected_genes and n_vars != expected_genes:
            issues.append(f"genes: expected {expected_genes}, got {n_vars}")
            ok = False

        return {
            "valid": ok,
            "n_cells": n_obs,
            "n_genes": n_vars,
            "issues": issues,
        }
    except Exception as e:
        return {"valid": False, "error": str(e)}


def validate_rds(filepath, expected_cells, expected_genes):
    """用 crosscell inspect 或 R 验证 RDS 文件"""
    # 先尝试 crosscell inspect
    try:
        result = subprocess.run(
            ["crosscell", "inspect", "--input", filepath],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode == 0:
            import re
            output = result.stdout
            cells_match = re.search(r"Cells: (\d+)", output)
            genes_match = re.search(r"Genes: (\d+)", output)
            if cells_match and genes_match:
                n_cells = int(cells_match.group(1))
                n_genes = int(genes_match.group(1))
                ok = True
                issues = []
                if expected_cells and n_cells != expected_cells:
                    issues.append(f"cells: expected {expected_cells}, got {n_cells}")
                    ok = False
                if expected_genes and n_genes != expected_genes:
                    issues.append(f"genes: expected {expected_genes}, got {n_genes}")
                    ok = False
                return {"valid": ok, "n_cells": n_cells, "n_genes": n_genes, "issues": issues}
    except Exception:
        pass

    # 回退到 R
    try:
        script = f"""
obj <- readRDS("{filepath}")
if (inherits(obj, "Seurat")) {{
    cat(ncol(obj), ",", nrow(obj), "\\n", sep="")
}} else {{
    cat("unknown_type\\n")
}}
"""
        result = subprocess.run(
            ["Rscript", "-e", script],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            return {"valid": False, "error": result.stderr[:300]}

        line = result.stdout.strip()
        if line == "unknown_type":
            return {"valid": False, "error": "not a Seurat object"}

        parts = line.split(",")
        n_cells, n_genes = int(parts[0]), int(parts[1])
        ok = True
        issues = []
        if expected_cells and n_cells != expected_cells:
            issues.append(f"cells: expected {expected_cells}, got {n_cells}")
            ok = False
        if expected_genes and n_genes != expected_genes:
            issues.append(f"genes: expected {expected_genes}, got {n_genes}")
            ok = False
        return {"valid": ok, "n_cells": n_cells, "n_genes": n_genes, "issues": issues}
    except Exception as e:
        return {"valid": False, "error": str(e)}


def build_expected_map(matrix):
    """从 test_matrix 构建 test_id → (expected_cells, expected_genes) 映射"""
    expected = {}
    for direction in ("rds_to_h5ad", "h5ad_to_rds"):
        for case in matrix["test_cases"].get(direction, []):
            expected[case["id"]] = {
                "n_cells": case.get("n_cells"),
                "n_genes": case.get("n_genes"),
                "direction": direction,
            }
    return expected


def validate_tool(tool_name, matrix):
    """验证一个工具的所有输出"""
    results_data = load_results(tool_name)
    if not results_data:
        print(f"  ⚠️  没有找到 {tool_name} 的结果文件")
        return None

    expected_map = build_expected_map(matrix)
    validation = {"tool": tool_name, "results": [], "summary": {}}
    total = 0
    valid_count = 0
    invalid_count = 0
    skipped = 0

    for direction in ("rds_to_h5ad", "h5ad_to_rds"):
        for entry in results_data.get(direction, []):
            test_id = entry["test_id"]
            if entry["status"] != "success":
                continue

            total += 1
            ext = "h5ad" if direction == "rds_to_h5ad" else "rds"
            output_path = os.path.join(OUTPUT_TMP_DIR, f"{tool_name}_{test_id}.{ext}")

            if not os.path.exists(output_path):
                skipped += 1
                validation["results"].append({
                    "test_id": test_id, "direction": direction,
                    "status": "missing", "reason": "output file not found",
                })
                continue

            exp = expected_map.get(test_id, {})
            exp_cells = exp.get("n_cells")
            exp_genes = exp.get("n_genes")

            if ext == "h5ad":
                v = validate_h5ad(output_path, exp_cells, exp_genes)
            else:
                v = validate_rds(output_path, exp_cells, exp_genes)

            if v.get("valid"):
                valid_count += 1
                status = "valid"
                print(f"    ✅ {test_id}: {v.get('n_cells')} × {v.get('n_genes')}")
            else:
                invalid_count += 1
                status = "invalid"
                issues = v.get("issues", [v.get("error", "unknown")])
                print(f"    ❌ {test_id}: {'; '.join(str(i) for i in issues)}")

            validation["results"].append({
                "test_id": test_id,
                "direction": direction,
                "status": status,
                **v,
            })

    validation["summary"] = {
        "total": total,
        "valid": valid_count,
        "invalid": invalid_count,
        "skipped": skipped,
    }
    return validation


def main():
    parser = argparse.ArgumentParser(description="Validate benchmark output files")
    parser.add_argument("--tool", required=True, help="Tool name or 'all'")
    parser.add_argument("--matrix", default=MATRIX_FILE, help="test_matrix.json path")
    parser.add_argument("--output-dir", default=RESULTS_DIR, help="Results directory")
    args = parser.parse_args()

    matrix = load_matrix(args.matrix)

    if args.tool == "all":
        tools = [t["name"] for t in matrix["tools"]]
    else:
        tools = [args.tool]

    all_validations = {}
    for tool in tools:
        print(f"\n{'='*50}")
        print(f"  验证: {tool}")
        print(f"{'='*50}")
        v = validate_tool(tool, matrix)
        if v:
            all_validations[tool] = v
            s = v["summary"]
            print(f"\n  📊 {tool}: {s['valid']}/{s['total']} 有效, "
                  f"{s['invalid']} 无效, {s['skipped']} 缺失")

    # 保存验证结果
    output_file = os.path.join(args.output_dir, "validation.json")
    with open(output_file, "w") as f:
        json.dump(all_validations, f, indent=2, ensure_ascii=False)
    print(f"\n💾 验证结果已保存: {output_file}")


if __name__ == "__main__":
    main()
