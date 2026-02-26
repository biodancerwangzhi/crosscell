#!/usr/bin/env python3
"""
Export benchmark JSON results to CSV for supplementary materials.

Usage:
    python benchmark/scripts/export_results_to_csv.py

Output directory: benchmark/results/csv/
"""

import json
import csv
import os
from pathlib import Path

RESULTS_DIR = Path("benchmark/results")
CSV_DIR = RESULTS_DIR / "csv"


def load_json(filename):
    path = RESULTS_DIR / filename
    if not path.exists():
        return None
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def write_csv(filename, rows, fieldnames):
    path = CSV_DIR / filename
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  ✅ {filename} ({len(rows)} rows)")


def export_robustness():
    """Supp Table 1: Robustness benchmark — all tools, all datasets."""
    data = load_json("robustness_benchmark.json")
    if not data:
        return
    rows = []
    for tool, tool_data in data.items():
        for direction in ["rds_to_h5ad", "h5ad_to_rds"]:
            for r in tool_data.get(direction, []):
                rows.append({
                    "tool": tool,
                    "direction": direction,
                    "test_id": r.get("test_id", ""),
                    "dataset": r.get("dataset", ""),
                    "file": r.get("file", ""),
                    "seurat_version": r.get("seurat_version", ""),
                    "type": r.get("type", ""),
                    "scale": r.get("scale", ""),
                    "n_cells": r.get("n_cells", ""),
                    "n_genes": r.get("n_genes", ""),
                    "spatial": r.get("spatial", False),
                    "status": r.get("status", ""),
                    "conversion_time_seconds": r.get("conversion_time_seconds", ""),
                    "peak_memory_mb": r.get("peak_memory_mb", ""),
                    "error": r.get("error_message", r.get("error", "")),
                })
    fields = ["tool", "direction", "test_id", "dataset", "file",
              "seurat_version", "type", "scale", "n_cells", "n_genes",
              "spatial", "status", "conversion_time_seconds", "peak_memory_mb", "error"]
    write_csv("supp_table1_robustness.csv", rows, fields)


def export_tool_results():
    """Per-tool detailed results (crosscell, zellkonverter, anndataR, convert2anndata, easySCF)."""
    tool_files = {
        "crosscell": "crosscell.json",
        "zellkonverter": "zellkonverter.json",
        "anndataR": "anndataR.json",
        "convert2anndata": "convert2anndata.json",
        "easySCF": "easySCF.json",
    }
    for tool_name, filename in tool_files.items():
        data = load_json(filename)
        if not data:
            continue
        rows = []
        for direction in ["rds_to_h5ad", "h5ad_to_rds"]:
            for r in data.get(direction, []):
                rows.append({
                    "tool": data.get("tool", tool_name),
                    "direction": direction,
                    "test_id": r.get("test_id", ""),
                    "dataset": r.get("dataset", ""),
                    "file": r.get("file", ""),
                    "seurat_version": r.get("seurat_version", ""),
                    "type": r.get("type", ""),
                    "n_cells": r.get("n_cells", ""),
                    "n_genes": r.get("n_genes", ""),
                    "status": r.get("status", ""),
                    "conversion_time_seconds": r.get("conversion_time_seconds", ""),
                    "peak_memory_mb": r.get("peak_memory_mb", ""),
                    "output_size_mb": r.get("output_size_mb", ""),
                    "error": r.get("error_message", r.get("error", "")),
                })
        fields = ["tool", "direction", "test_id", "dataset", "file",
                  "seurat_version", "type", "n_cells", "n_genes",
                  "status", "conversion_time_seconds", "peak_memory_mb",
                  "output_size_mb", "error"]
        write_csv(f"tool_{tool_name}.csv", rows, fields)


def export_roundtrip_fidelity():
    """Roundtrip fidelity results."""
    data = load_json("roundtrip_fidelity.json")
    if not data:
        return

    # RDS roundtrip
    rows = []
    for r in data.get("rds_roundtrip", []):
        rows.append({
            "test_id": r.get("test_id", ""),
            "dataset": r.get("dataset", ""),
            "n_cells": r.get("n_cells", ""),
            "n_genes": r.get("n_genes", ""),
            "status": r.get("status", ""),
            "pearson_r": r.get("pearson_r", ""),
            "mse": r.get("mse", ""),
            "exact_match": r.get("exact_match", ""),
            "obs_match": r.get("obs_match", ""),
            "var_match": r.get("var_match", ""),
        })
    if rows:
        fields = ["test_id", "dataset", "n_cells", "n_genes", "status",
                  "pearson_r", "mse", "exact_match", "obs_match", "var_match"]
        write_csv("roundtrip_rds.csv", rows, fields)

    # H5AD roundtrip
    rows = []
    for r in data.get("h5ad_roundtrip", []):
        rows.append({
            "test_id": r.get("test_id", ""),
            "dataset": r.get("dataset", ""),
            "n_cells": r.get("n_cells", ""),
            "n_genes": r.get("n_genes", ""),
            "status": r.get("status", ""),
            "pearson_r": r.get("pearson_r", ""),
            "mse": r.get("mse", ""),
            "exact_match": r.get("exact_match", ""),
            "obs_match": r.get("obs_match", ""),
            "var_match": r.get("var_match", ""),
        })
    if rows:
        write_csv("roundtrip_h5ad.csv", rows, fields)

    # Mixed roundtrip
    mixed = data.get("mixed_roundtrip", [])
    rows = []
    if isinstance(mixed, list):
        for r in mixed:
            rows.append({
                "tool": r.get("tool", ""),
                "test_id": r.get("test_id", ""),
                "dataset": r.get("dataset", ""),
                "status": r.get("status", ""),
                "pearson_r": r.get("pearson_r", ""),
                "mse": r.get("mse", ""),
                "exact_match_pct": r.get("exact_match_pct", r.get("exact_match", "")),
                "error": r.get("error_message", r.get("error", "")),
            })
    elif isinstance(mixed, dict):
        for tool, results in mixed.items():
            for r in results:
                rows.append({
                    "tool": tool,
                    "test_id": r.get("test_id", ""),
                    "dataset": r.get("dataset", ""),
                    "status": r.get("status", ""),
                    "pearson_r": r.get("pearson_r", ""),
                    "mse": r.get("mse", ""),
                    "exact_match_pct": r.get("exact_match_pct", r.get("exact_match", "")),
                    "error": r.get("error_message", r.get("error", "")),
                })
    if rows:
        fields_mixed = ["tool", "test_id", "dataset", "status",
                        "pearson_r", "mse", "exact_match_pct", "error"]
        write_csv("roundtrip_mixed.csv", rows, fields_mixed)


def export_type_fidelity():
    """Type fidelity results."""
    data = load_json("type_fidelity.json")
    if not data:
        return

    type_results = data.get("type_results", {})
    expected = data.get("expected_types", {})

    if isinstance(type_results, dict):
        rows = []
        for tool, results in type_results.items():
            if isinstance(results, dict):
                for col, actual_type in results.items():
                    rows.append({
                        "tool": tool,
                        "column": col,
                        "expected_type": expected.get(col, ""),
                        "actual_type": actual_type,
                        "preserved": actual_type == expected.get(col, ""),
                    })
        if rows:
            fields = ["tool", "column", "expected_type", "actual_type", "preserved"]
            write_csv("type_fidelity.csv", rows, fields)

    # Retention summary
    retention = data.get("retention", {})
    if retention:
        rows2 = [{"tool": t, "retention_pct": v} for t, v in retention.items()]
        write_csv("type_retention_summary.csv", rows2, ["tool", "retention_pct"])


def export_multi_tool_fidelity():
    """Multi-tool fidelity comparison."""
    data = load_json("multi_tool_fidelity.json")
    if not data:
        return

    # CrossCell roundtrip
    rows = []
    for r in data.get("crosscell_roundtrip", []):
        rows.append({
            "direction": r.get("direction", ""),
            "test_id": r.get("test_id", ""),
            "dataset": r.get("dataset", ""),
            "status": r.get("status", ""),
            "pearson_r": r.get("pearson_r", ""),
            "mse": r.get("mse", ""),
            "exact_match": r.get("exact_match", ""),
            "max_abs_diff": r.get("max_abs_diff", ""),
        })
    if rows:
        fields = ["direction", "test_id", "dataset", "status",
                  "pearson_r", "mse", "exact_match", "max_abs_diff"]
        write_csv("fidelity_crosscell_roundtrip.csv", rows, fields)

    # Cross-tool comparison
    rows = []
    for r in data.get("cross_tool", []):
        rows.append({
            "tool": r.get("tool", ""),
            "test_id": r.get("test_id", ""),
            "dataset": r.get("dataset", ""),
            "status": r.get("status", ""),
            "pearson_r": r.get("pearson_r", ""),
            "mse": r.get("mse", ""),
            "exact_match": r.get("exact_match", ""),
            "max_abs_diff": r.get("max_abs_diff", ""),
        })
    if rows:
        fields2 = ["tool", "test_id", "dataset", "status",
                   "pearson_r", "mse", "exact_match", "max_abs_diff"]
        write_csv("fidelity_cross_tool.csv", rows, fields2)


def export_cellxgene_scalability():
    """CELLxGENE scalability results."""
    data = load_json("cellxgene_scalability.json")
    if not data:
        return
    rows = []
    for r in data.get("results", []):
        rows.append({
            "dataset": r.get("dataset", ""),
            "n_cells": r.get("n_cells", ""),
            "n_genes": r.get("n_genes", ""),
            "file_size_gb": r.get("file_size_gb", r.get("file_size_mb", "")),
            "inspect_status": r.get("inspect_status", r.get("inspect", "")),
            "h2r_status": r.get("h2r_status", r.get("h5ad_to_rds", "")),
            "h2r_time": r.get("h2r_time", ""),
            "r2h_status": r.get("r2h_status", r.get("rds_to_h5ad", "")),
            "r2h_time": r.get("r2h_time", ""),
            "peak_memory_gb": r.get("peak_memory_gb", ""),
            "test_env": r.get("test_env", r.get("environment", "")),
        })
    if rows:
        fields = ["dataset", "n_cells", "n_genes", "file_size_gb",
                  "inspect_status", "h2r_status", "h2r_time",
                  "r2h_status", "r2h_time", "peak_memory_gb", "test_env"]
        write_csv("cellxgene_scalability.csv", rows, fields)


def export_tool_versions():
    """Tool versions."""
    data = load_json("tool_versions.json")
    if not data:
        return
    tools = data.get("tools", {})
    rows = []
    if isinstance(tools, dict):
        for name, info in tools.items():
            rows.append({
                "tool": name,
                "version": info.get("version", ""),
                "language": info.get("language", ""),
                "type": info.get("type", ""),
                "note": info.get("note", ""),
            })
    elif isinstance(tools, list):
        for tool in tools:
            rows.append({
                "tool": tool.get("name", ""),
                "version": tool.get("version", ""),
                "language": tool.get("language", ""),
                "type": tool.get("type", ""),
                "note": tool.get("note", ""),
            })
    if rows:
        write_csv("tool_versions.csv", rows, ["tool", "version", "language", "type", "note"])


def export_dataset_summary():
    """Dataset summary."""
    data = load_json("dataset_summary.json")
    if not data:
        return
    rows = []
    for d in data.get("datasets", []):
        rows.append({
            "name": d.get("name", ""),
            "format": d.get("format", ""),
            "n_cells": d.get("n_cells", ""),
            "n_genes": d.get("n_genes", ""),
            "scale": d.get("scale", ""),
            "source": d.get("source", ""),
        })
    if rows:
        write_csv("dataset_summary.csv", rows,
                  ["name", "format", "n_cells", "n_genes", "scale", "source"])


def export_tools_fidelity():
    """Tools fidelity summary (from tools_fidelity.json)."""
    data = load_json("tools_fidelity.json")
    if not data or not isinstance(data, list):
        return
    if len(data) == 0:
        return
    # Infer fields from first record
    fields = list(data[0].keys())
    write_csv("tools_fidelity_summary.csv", data, fields)


def export_memory_scaling():
    """Memory scaling results."""
    data = load_json("memory_scaling.json")
    if not data:
        return
    rows = []
    for r in data.get("results", []):
        row = {
            "n_cells": r.get("n_cells", ""),
            "nnz": r.get("nnz", ""),
            "peak_memory_mb": r.get("peak_memory_mb", ""),
            "time_seconds": r.get("time_seconds", ""),
        }
        rows.append(row)
    if rows:
        fields = ["n_cells", "nnz", "peak_memory_mb", "time_seconds"]
        write_csv("memory_scaling.csv", rows, fields)


def export_h5ad_read_tests():
    """anndataR and easySCF H5AD read test results."""
    # anndataR
    data = load_json("anndataR_h5ad_read_test.json")
    if data and isinstance(data, dict):
        rows = []
        for dataset, result in data.items():
            if isinstance(result, dict):
                rows.append({
                    "dataset": dataset,
                    "read_h5ad": result.get("read_h5ad", ""),
                    "to_Seurat": result.get("to_Seurat", ""),
                    "error": result.get("error", ""),
                })
        if rows:
            write_csv("anndataR_h5ad_read_test.csv", rows,
                      ["dataset", "read_h5ad", "to_Seurat", "error"])

    # easySCF
    data = load_json("easyscf_h5ad_to_rds_test.json")
    if data and isinstance(data, list):
        if len(data) > 0:
            fields = list(data[0].keys())
            write_csv("easyscf_h5ad_to_rds_test.csv", data, fields)


def export_nnz():
    """NNZ (non-zero elements) data for H5AD and RDS files."""
    for name in ["nnz_h5ad", "nnz_rds"]:
        data = load_json(f"{name}.json")
        if not data or not isinstance(data, dict):
            continue
        rows = []
        for filename, info in data.items():
            if isinstance(info, dict):
                row = {"file": filename}
                row.update(info)
                rows.append(row)
            else:
                rows.append({"file": filename, "nnz": info})
        if rows:
            fields = list(rows[0].keys())
            write_csv(f"{name}.csv", rows, fields)


def main():
    os.makedirs(CSV_DIR, exist_ok=True)
    print(f"Exporting benchmark results to {CSV_DIR}/\n")

    print("1. Robustness benchmark (Supp Table 1)")
    export_robustness()

    print("\n2. Per-tool detailed results")
    export_tool_results()

    print("\n3. Roundtrip fidelity")
    export_roundtrip_fidelity()

    print("\n4. Type fidelity")
    export_type_fidelity()

    print("\n5. Multi-tool fidelity comparison")
    export_multi_tool_fidelity()

    print("\n6. CELLxGENE scalability")
    export_cellxgene_scalability()

    print("\n7. Tool versions")
    export_tool_versions()

    print("\n8. Dataset summary")
    export_dataset_summary()

    print("\n9. Tools fidelity summary")
    export_tools_fidelity()

    print("\n10. Memory scaling")
    export_memory_scaling()

    print("\n11. H5AD read tests (anndataR, easySCF)")
    export_h5ad_read_tests()

    print("\n12. NNZ data")
    export_nnz()

    print(f"\n✅ All CSV files exported to {CSV_DIR}/")


if __name__ == "__main__":
    main()
