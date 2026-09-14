# R 环境搭建记录（第一步完成）

日期：2026-09-13　｜　机器：Windows 10 (19045)，非管理员账户，8 核

---

## 1. 结论

**R 4.6.1 已在本机可用，scTenifoldKnk 全流程已跑通。** 你可以直接开始用，不需要再做任何安装。

| 项目 | 结果 |
|---|---|
| R 版本 | 4.6.1 (2026-06-24 ucrt) |
| 安装位置 | `work/R/`（项目目录内，可整体移动或删除） |
| R 库路径 | `work/R/library` |
| 已安装包 | **49 个**（含 scTenifoldKnk 1.1、scTenifoldNet 1.4、locfdr 1.1-8） |
| 磁盘占用 | work 目录合计约 519 MB（R 323 MB + 依赖包 74 MB + 安装器 107 MB） |
| 冒烟测试 | 通过（100 基因 × 500 细胞，0.11 分钟） |
| 面板脚本自测 | 通过（20 基因面板 + 稳定性 + 负控 + 出图） |

---

## 2. 为什么不是"双击安装"

标准安装路径全部被系统权限挡住，逐一记录如下：

| 尝试 | 结果 |
|---|---|
| `winget install RProject.R` | winget.exe 在沙箱内无法启动（WindowsApps 应用别名被拒） |
| 官方安装器静默安装（`/VERYSILENT /DIR=...`） | **ExitCode 5**：文件写入成功，但创建开始菜单快捷方式与写入注册表 `HKCU\Software\R-core` 被拒（错误 5），安装器自动回滚 |
| 最终采用：**解包式安装** | 用 `innoextract 1.9` 直接解出官方安装器内的全部文件，绕过注册表与快捷方式 |

解包安装得到的是**完整、原版的 CRAN R**（不是精简版），代价是：

- 没有开始菜单/桌面快捷方式（用命令行启动即可）
- 没有注册表项（R 依据 `Rscript.exe` 的路径自动推断 `R_HOME`，实测正常）
- 未加入系统 PATH（用完整路径调用即可）

---

## 3. 怎么用

所有命令都在项目根目录 `C:\Users\realb\Documents\Codex\2026-09-13\sctenifoldknk-knockdown-x20` 下执行。

```powershell
# 版本确认
.\work\R\bin\Rscript.exe -e "cat(R.version.string)"

# 环境体检 + 冒烟测试（离线模式会自动跳过联网安装）
.\work\R\bin\Rscript.exe outputs\scripts\00_setup_environment.R

# 面板虚拟敲除主脚本（参数可在 work/vko_overrides.R 里覆盖）
.\work\R\bin\Rscript.exe outputs\scripts\run_vko_panel.R
```

开交互式 R 会话：

```powershell
.\work\R\bin\R.exe
```

> 想把 R 移到别处（例如 `C:\Users\realb\R\R-4.6.1`）直接整体移动 `work/R` 即可，
> R 会自动重新定位 `R_HOME`，只需相应修改上面的命令路径。

---

## 4. 本机网络限制（重要）

**R 在本机无法联网**，已实测确认：

```
download.file(method = "auto")    -> FAIL: SSL connect error
download.file(method = "libcurl") -> FAIL: SSL connect error (schannel: SEC_E_NO_CREDENTIALS)
download.file(method = "wininet") -> FAIL: InternetOpenUrl failed
download.file(method = "curl")    -> FAIL: curl exit 35, schannel 凭证错误
```

原因是沙箱环境下 Windows 的 schannel 证书凭据不可用；Node.js 自带 OpenSSL，不受影响。

因此本环境的 R 包安装采用**离线两段式**：

1. 用 Node 抓依赖树和二进制包（Node 能联网）
   `node work\fetch_r_pkgs.mjs`
2. 用 R 从本地 zip 离线安装
   `.\work\R\bin\Rscript.exe work\install_r_pkgs.R`

> 影响：`plotKO(annotate = TRUE)` 与 `enrichR` 相关功能需要联网访问 Enrichr，本机用不了。
> 主流程（建网、虚拟敲除、差异调控、稳定性、阴性对照）**完全不受影响**。

---

## 5. 已安装的 49 个包

核心：`scTenifoldKnk 1.1`、`scTenifoldNet 1.4`、`Matrix 1.7.5`、`MASS 7.3-65`、`locfdr 1.1-8`

依赖链：`Rcpp 1.1.2`、`RcppEigen 0.3.4.0.2`、`RcppArmadillo 15.6.0-1`、`RSpectra 0.16-2`、`RhpcBLASctl`、`cli 3.6.6`、`enrichR 3.4`、`igraph 2.3.3`、`reshape2 1.4.5`、`plyr`、`stringi`、`stringr`、`ggplot2 4.0.3`、`future`/`furrr`/`globals`/`listenv`/`parallelly`、`httr`/`curl`/`openssl`/`jsonlite`/`rjson`/`askpass`/`sys`/`mime`、`rlang`/`vctrs`/`glue`/`magrittr`/`lifecycle`/`purrr`/`withr`/`R6`/`S7`、`scales`/`gtable`/`isoband`/`farver`/`labeling`/`RColorBrewer`/`viridisLite`、`WriteXLS`、`cpp11`、`digest`、`pkgconfig`

**未安装**（属 Bioconductor 或数据准备阶段才需要）：

| 包 | 来源 | 影响 |
|---|---|---|
| `fgsea` | Bioconductor | 富集分析会被自动跳过（脚本已优雅降级） |
| `preprocessCore` | Bioconductor | 稳定性分析的分位数归一化被跳过（已优雅降级） |
| `Seurat` / `SeuratObject` | CRAN | 数据准备阶段（W1 注释脚本）需要 |

---

## 6. 自测记录

### 6.1 官方示例冒烟测试

```
输入：100 基因 × 500 细胞（负二项模拟，QC 后 249 细胞）
参数：nc_nNet = 5, nc_nCells = 300, td_K = 3
结果：用时 0.11 分钟；输出列完整（gene/distance/Z/FC/p.value/p.adj）
      Top1 = ng10（被敲除的基因本身），FDR = 3.0e-20   ← 合理性检查通过
```

### 6.2 面板脚本端到端自测

用 `work/make_test_data.R` 生成合成数据（1200 基因 × 900 细胞，含 Microglia/Astrocyte 两种细胞类型、Control/Disease 两种条件）验证 `run_vko_panel.R`：

- 数据读取、细胞筛选、QC、核糖体/线粒体过滤 ✅
- `transcriptomeWide` 面板模式（一次建网 + 20 个基因扰动）✅
- 距离 → Z / FC / p / FDR 换算 ✅
- 面板汇总、Top 命中表、稳定性分析、阴性对照比较 ✅
- 单基因模式 + `plotKO()` 子网络出图 ✅
- 产物：距离矩阵、全量统计表、Top 命中表、面板汇总、稳定性表、子网络 PNG、参数与 sessionInfo ✅

产出目录：`work/test_output/`（自测用，非正式结果）

> 自测用的是随机合成数据，几乎没有真实调控结构，所以每个基因只有 1 个显著扰动基因（即被敲除基因本身），符合预期。真实数据上会看到成百上千个被扰动基因。

---

## 7. 在你自己的机器或实验室服务器上安装

若那台机器有管理员权限且能联网，标准安装即可，**不需要**上面这套解包流程：

```r
install.packages("scTenifoldKnk")           # CRAN v1.1
install.packages(c("locfdr", "igraph", "reshape2", "enrichR"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("fgsea", "preprocessCore"))
```

然后把 `outputs/scripts/` 下的脚本拷过去即可，注意修改配置区的路径。

若那台机器同样无管理员权限或无法联网，可复刻本机流程：
把 `work/installers/R-4.6.1-win.exe`、`work/installers/innoextract/`、`work/fetch_r_pkgs.mjs`、`work/install_r_pkgs.R` 一起拷过去，按第 2、4 节操作。

---

## 8. 第一步完成状态

| 检查点 | 状态 |
|---|---|
| R 可运行 | ✅ 4.6.1 |
| scTenifoldKnk ≥ 1.1 可载入 | ✅ 1.1，4 个导出函数齐全 |
| 官方示例跑通 | ✅ |
| 面板脚本端到端跑通 | ✅ |
| 掌握本机耗时基线 | ✅ 见 `outputs/R_benchmark_CN.md` |
| CRAN 直连 | ❌ 本机不可用（已用离线流程绕过） |

**下一步（W1）：选数据。** 打开 `outputs/templates/dataset_survey.csv` 选定主数据集与复现数据集，
再按 `work/make_test_data.R` 的格式把真实数据整理成 `work/data/` 下的四个文件即可直接开跑。
