# CHARM Preprocess Pipeline

**概述**
- 本流程用于 CHARM 原始数据预处理，包含 RNA、Hi-C 和 CUT&Tag/ATAC 等模块。

**目录结构**
- `work_dir` 是数据集根目录，必须包含 `Rawdata/`。
- 期望的输入结构：`Rawdata/<sample>/<sample>_1.fq.gz` 和 `Rawdata/<sample>/<sample>_2.fq.gz`。
- 运行后会在 `work_dir` 下生成 `processed/`、`result/`、`stat/`、`analysis/`。

**配置**
- 配置文件：`CHARM_preprocess_pipeline/config.yaml`。
- `work_dir` 可写为绝对或相对路径，相对路径会以流程目录为基准解析。
- 资源参数统一在 `resources` 中设置（仅 `threads`/`mem_mb`/`gpu_proc`）。
- 通过 `if_*` 和 `ref_genome` 控制模块开关与参考版本。

**运行**
- 推荐：`bash CHARM_preprocess_pipeline/runCHARM.sh`
- 或手动：`snakemake -j 190 --resources gpu_proc=20 -s CHARM_preprocess_pipeline/CHARM.smk --rerun-incomplete --keep-going `

**输出**
- `processed/` 中间文件与对齐结果。
- `result/` 主要结果文件。
- `stat/` 统计信息。
- `analysis/` 自动拷贝的 `stat.ipynb`。
# CHARM_preprocess_pipeline
