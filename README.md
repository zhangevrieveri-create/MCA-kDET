# MCA-kDET 模拟与实证代码

- `simulation/`：主模拟。`run.R`（入口，`Rscript run.R --exp=E1 --profile=quick --cores=30`）、`R/`（模型、学习器、置换引擎、实验网格；`models.R`、`calib_params*.R`、`tests.R` 已冻结）、`report_*.R`（出表）、`tests/`（单元测试与诊断）、`run_*.bat`（一键运行）。
- `extra_sims/`：第三轮功效/形状补充模拟（motivated_v4、靶向模拟），仅脚本与设计文件。
- `applications/`：单细胞与金融实证脚本（原始数据、结果未入库）。
- `figure_table_scripts/`：生成论文图表和补充表的脚本。
- `docs/`：交接、设计、运行说明与报告。

在 `simulation/` 下运行；R 路径为 `D:/R-4.6.1/bin/Rscript.exe`（不稳时用 4.3.3）。
