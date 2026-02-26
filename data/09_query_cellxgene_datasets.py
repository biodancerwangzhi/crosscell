#!/usr/bin/env python3
"""
CELLxGENE Benchmark 数据集选择文档

本脚本记录了 10 个 CELLxGENE benchmark 数据集的选择标准和来源信息。
这些数据集用于 CrossCell 论文的大规模 scalability 测试（Fig 5）。

选择策略：分层采样（stratified sampling）
  1. 组织多样性：10 种组织（blood, lung, heart, kidney, gut, liver, brain×3, pancreas, skin）
  2. 技术多样性：10x 3'/5' v2/v3, 10x multiome, scATAC-seq, Smart-seq2
  3. 规模梯度：10k → 122k → 428k → 836k cells
  4. 疾病状态：normal + disease (COVID-19, BCC tumor)

排除说明：
  - retina_244k：18 GB 文件含异常 dense 矩阵结构，inspect 需 62 GB 内存，非工具问题

数据来源：
  - 所有数据集均来自 CELLxGENE Discover (https://cellxgene.cziscience.com)
  - 每个数据集有永久 URL，格式为 https://datasets.cellxgene.cziscience.com/{dataset_id}.h5ad
  - 初始候选列表通过 CELLxGENE 网页浏览器筛选（Homo sapiens, >10k cells）

下载方式：
  python3 data/05_download_h5ad_datasets.py --source cellxgene-benchmark

用法（仅输出数据集信息，不执行下载）：
  python3 data/09_query_cellxgene_datasets.py
"""

import json

# 10 个最终选定的 CELLxGENE benchmark 数据集
CELLXGENE_BENCHMARK_DATASETS = [
    {
        "name": "cellxgene_skin_bcc_10k",
        "cells": 9841,
        "genes": 26886,
        "tissue": "skin",
        "assay": "multi-assay",
        "disease": "basal cell carcinoma",
        "sparsity": "93.1%",
        "file_size_mb": 89,
        "scale": "small",
        "dataset_id": "a6670528-1e94-4920-a78f-49a76b96ca9e",
        "url": "https://datasets.cellxgene.cziscience.com/a6670528-1e94-4920-a78f-49a76b96ca9e.h5ad",
        "selection_reason": "Small tumor dataset, disease state diversity",
    },
    {
        "name": "cellxgene_tabula_liver_22k",
        "cells": 22214,
        "genes": 60606,
        "tissue": "liver",
        "assay": "Smart-seq2 + 10x 3' v3",
        "disease": "normal",
        "sparsity": "95.4%",
        "file_size_mb": 796,
        "scale": "medium",
        "dataset_id": "62e4e805-9a9a-42ca-a6e7-85df258b382c",
        "url": "https://datasets.cellxgene.cziscience.com/62e4e805-9a9a-42ca-a6e7-85df258b382c.h5ad",
        "selection_reason": "Tabula Sapiens reference, Smart-seq2 technology diversity",
    },
    {
        "name": "cellxgene_kidney_atacseq_37k",
        "cells": 37747,
        "genes": 19276,
        "tissue": "kidney",
        "assay": "10x scATAC-seq",
        "disease": "normal",
        "sparsity": "69.4%",
        "file_size_mb": 898,
        "scale": "medium",
        "dataset_id": "fa0f9f5e-e7fb-42a4-a801-3410b0923ccc",
        "url": "https://datasets.cellxgene.cziscience.com/fa0f9f5e-e7fb-42a4-a801-3410b0923ccc.h5ad",
        "selection_reason": "scATAC-seq technology, low sparsity (69.4%) stress test",
    },
    {
        "name": "cellxgene_brain_multiome_102k",
        "cells": 101924,
        "genes": 35451,
        "tissue": "brain",
        "assay": "10x multiome",
        "disease": "normal",
        "sparsity": "93.2%",
        "file_size_mb": 1142,
        "scale": "large",
        "dataset_id": "412352dd-a919-4d8e-9f74-e210627328b5",
        "url": "https://datasets.cellxgene.cziscience.com/412352dd-a919-4d8e-9f74-e210627328b5.h5ad",
        "selection_reason": "10x multiome technology, 100k+ scale threshold",
    },
    {
        "name": "cellxgene_pancreas_122k",
        "cells": 121916,
        "genes": 32356,
        "tissue": "pancreas",
        "assay": "multi-assay",
        "disease": "normal",
        "sparsity": "96.0%",
        "file_size_mb": 915,
        "scale": "large",
        "dataset_id": "72e6a74d-8804-4305-81d7-a12b93f2ac0d",
        "url": "https://datasets.cellxgene.cziscience.com/72e6a74d-8804-4305-81d7-a12b93f2ac0d.h5ad",
        "selection_reason": "HCA harmonized atlas, 16 GB machine upper limit (~122k cells)",
    },
    {
        "name": "cellxgene_brain_dlpfc_172k",
        "cells": 172120,
        "genes": 37490,
        "tissue": "brain",
        "assay": "10x 3' v3",
        "disease": "normal",
        "sparsity": "91.9%",
        "file_size_mb": 2234,
        "scale": "very_large",
        "dataset_id": "9198437d-f781-4215-bf4f-32c3177e57df",
        "url": "https://datasets.cellxgene.cziscience.com/9198437d-f781-4215-bf4f-32c3177e57df.h5ad",
        "selection_reason": "First dataset requiring >16 GB RAM, 64 GB server test",
    },
    {
        "name": "cellxgene_gut_428k",
        "cells": 428469,
        "genes": 32383,
        "tissue": "gut",
        "assay": "10x 3'/5' v2",
        "disease": "normal",
        "sparsity": "94.6%",
        "file_size_mb": 5462,
        "scale": "very_large",
        "dataset_id": "f34d2b82-9265-4a73-bda4-852933bf2a8d",
        "url": "https://datasets.cellxgene.cziscience.com/f34d2b82-9265-4a73-bda4-852933bf2a8d.h5ad",
        "selection_reason": "Large intestinal atlas, 400k+ scale",
    },
    {
        "name": "cellxgene_heart_486k",
        "cells": 486134,
        "genes": 32383,
        "tissue": "heart",
        "assay": "10x 3' v2/v3",
        "disease": "normal",
        "sparsity": "96.1%",
        "file_size_mb": 2816,
        "scale": "very_large",
        "dataset_id": "a7f6822d-0e0e-451e-9858-af81392fcb9b",
        "url": "https://datasets.cellxgene.cziscience.com/a7f6822d-0e0e-451e-9858-af81392fcb9b.h5ad",
        "selection_reason": "Large cardiac atlas, near 500k scale",
    },
    {
        "name": "cellxgene_hlca_core_585k",
        "cells": 584944,
        "genes": 27402,
        "tissue": "lung",
        "assay": "multi-assay",
        "disease": "normal",
        "sparsity": None,
        "file_size_mb": 5602,
        "scale": "very_large",
        "dataset_id": "4cb45d80-499a-48ae-a056-c71ac3552c94",
        "url": "https://datasets.cellxgene.cziscience.com/4cb45d80-499a-48ae-a056-c71ac3552c94.h5ad",
        "selection_reason": "Human Lung Cell Atlas core, 500k+ scale",
    },
    {
        "name": "cellxgene_combat_pbmc_836k",
        "cells": 836148,
        "genes": 36306,
        "tissue": "blood",
        "assay": "10x 5' v1",
        "disease": "COVID-19/sepsis/flu",
        "sparsity": None,
        "file_size_mb": 5337,
        "scale": "very_large",
        "dataset_id": "687c09ff-731a-4e3d-ac07-4c29c33a6338",
        "url": "https://datasets.cellxgene.cziscience.com/687c09ff-731a-4e3d-ac07-4c29c33a6338.h5ad",
        "selection_reason": "Largest dataset, 800k+ cells, disease cohort",
    },
]

# 百万级数据集（scalability 极限测试）
CELLXGENE_MILLION_DATASETS = [
    {
        "name": "cellxgene_eqtl_autoimmune_1.2M",
        "cells": 1248980,
        "genes": None,
        "tissue": "blood",
        "assay": "10x 3' v2",
        "disease": "normal",
        "sparsity": None,
        "file_size_mb": 4000,
        "scale": "million",
        "dataset_id": "a3f5651f-cd1a-4d26-8165-74964b79b4f2",
        "url": "https://datasets.cellxgene.cziscience.com/a3f5651f-cd1a-4d26-8165-74964b79b4f2.h5ad",
        "selection_reason": "1.2M cells, smaller file size (~4GB) may fit in 64GB RAM for full conversion",
        "reference": "Single-cell eQTL mapping identifies cell type specific genetic control of autoimmune disease",
    },
    {
        "name": "cellxgene_immune_1.8M",
        "cells": 1821725,
        "genes": 32357,
        "tissue": "blood",
        "assay": "10x 3' v2",
        "disease": "normal",
        "sparsity": None,
        "file_size_mb": 9460,
        "scale": "million",
        "dataset_id": "39832a61-3df1-4f1f-971d-3c2a74539047",
        "url": "https://datasets.cellxgene.cziscience.com/39832a61-3df1-4f1f-971d-3c2a74539047.h5ad",
        "selection_reason": "1.8M cells, OOM on 64GB — demonstrates dgCMatrix serialization limit",
        "reference": "Immune cell census",
    },
]

# 排除的数据集
EXCLUDED_DATASETS = [
    {
        "name": "cellxgene_retina_244k",
        "cells": 244474,
        "genes": 30172,
        "tissue": "retina",
        "file_size_mb": 17986,
        "dataset_id": "d4f1a799-9616-4089-90e9-f3fe843f01e5",
        "url": "https://datasets.cellxgene.cziscience.com/d4f1a799-9616-4089-90e9-f3fe843f01e5.h5ad",
        "exclusion_reason": "Abnormal dense matrix structure. 18 GB file (vs ~2 GB expected for 244k cells). "
                           "Inspect alone requires 62 GB RAM. Not a CrossCell issue — the file stores "
                           "expression data as dense array instead of sparse, making it an outlier.",
    },
]


def main():
    print("=" * 70)
    print("CELLxGENE Benchmark Dataset Selection")
    print("=" * 70)
    
    print(f"\n📊 Selected: {len(CELLXGENE_BENCHMARK_DATASETS)} datasets")
    print(f"🚀 Million-scale: {len(CELLXGENE_MILLION_DATASETS)} datasets")
    print(f"❌ Excluded: {len(EXCLUDED_DATASETS)} datasets")
    
    # Summary by scale
    scales = {}
    for d in CELLXGENE_BENCHMARK_DATASETS:
        s = d["scale"]
        scales[s] = scales.get(s, 0) + 1
    print(f"\n📏 Scale distribution: {scales}")
    
    # Summary by tissue
    tissues = set(d["tissue"] for d in CELLXGENE_BENCHMARK_DATASETS)
    print(f"🧬 Tissues: {len(tissues)} ({', '.join(sorted(tissues))})")
    
    # Summary by assay
    assays = set(d["assay"] for d in CELLXGENE_BENCHMARK_DATASETS)
    print(f"🔬 Assays: {len(assays)} ({', '.join(sorted(assays))})")
    
    # Print table
    print(f"\n{'Name':<35} {'Cells':>10} {'Genes':>8} {'Size MB':>8} {'Tissue':<10} {'Scale':<12}")
    print("-" * 95)
    for d in CELLXGENE_BENCHMARK_DATASETS:
        print(f"{d['name']:<35} {d['cells']:>10,} {d['genes']:>8,} {d['file_size_mb']:>8,} {d['tissue']:<10} {d['scale']:<12}")
    
    print(f"\n❌ Excluded:")
    for d in EXCLUDED_DATASETS:
        print(f"  {d['name']}: {d['exclusion_reason']}")
    
    print(f"\n🚀 Million-scale datasets:")
    for d in CELLXGENE_MILLION_DATASETS:
        print(f"  {d['name']:<40} {d['cells']:>12,} cells  {d['file_size_mb']:>6,} MB  {d.get('reference', '')}")
    
    # Output as JSON
    output = {
        "selection_criteria": {
            "tissue_diversity": "10 tissues (blood, brain, gut, heart, kidney, liver, lung, pancreas, skin)",
            "technology_diversity": "10x 3'/5' v2/v3, 10x multiome, scATAC-seq, Smart-seq2",
            "scale_gradient": "10k → 122k → 428k → 836k → 1.2M → 1.8M cells",
            "disease_states": "normal + disease (COVID-19, BCC tumor, autoimmune)",
        },
        "datasets": CELLXGENE_BENCHMARK_DATASETS,
        "million_scale": CELLXGENE_MILLION_DATASETS,
        "excluded": EXCLUDED_DATASETS,
    }
    
    print(f"\n💾 JSON output:")
    print(json.dumps(output, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
