# 发布清单：GitHub 仓库 + Zenodo DOI

更新日期：2026-09-14 ｜ 对应 `docs/09_methods_and_limitations_draft_EN.md` 的
`Data and code availability` 与根目录 `CITATION.cff`

---

## 0. 四个值的落实情况（2026-09-14 已填）

| 信息 | 值 | 状态 |
|---|---|---|
| 作者 + ORCID | 朱子龙 / **Zhu Zilong**，ORCID `0000-0002-6955-0903` | ✅ 已写入 `CITATION.cff` |
| GitHub 仓库地址 | `https://github.com/zhuzilong1976/zhuzilong` | ⚠️ 按"账号 zhuzilong1976 / 仓库 zhuzilong"推断填入，**若用户名不同请改这一行** |
| 论文题目 | `Containment, not biology: the differentially regulated gene set of single-cell virtual knockout is the knockout's own out-neighbourhood` | ✅ 已填入 `preferred-citation`（暂定） |
| 目标期刊 | `PLOS Computational Biology` | ⚠️ 暂定值，投稿前确认（见第 7 节的期刊建议） |
| git 提交身份 | `zhuzilong` / `zhuzilong1976@gmail.com` | ✅ 已配置为仓库级 `user.name`/`user.email` 并完成首次提交 |

只剩 Zenodo DOI 一项（第 3 节）。

---

## 1. 本地仓库状态

仓库已在 `vko-project/` 下初始化为 git 仓库，并已完成首次提交：

```text
37b0cc3  Initial release: code, figures and tables for the virtual-knockout containment study
         109 files, 0 files > 1 MB, working tree clean
```

尚未 `push`（需要你的 GitHub 凭据）。

`.gitignore` 已排除大文件与派生数据：`work/`、`outputs/`、`data/*`、`*.rds`、
`*.h5ad`、`*.pkl`、`*.mtx`、`*.zip`。因此入 git 的只有代码、文档、图与表
（约 1 MB）；完整距离矩阵等大文件走 Zenodo（第 3 节）。

## 2. 建立 GitHub 仓库并推送

```bash
cd C:\Users\朱子龙\Desktop\codex\vko-project

git branch -M main

# 先在 GitHub 网页新建一个空仓库（不要勾选 README/.gitignore），然后：
git remote add origin https://github.com/zhuzilong1976/zhuzilong.git
git push -u origin main
```

推送前用 `git status` 确认没有 `work/`、`outputs/` 下的文件。

## 3. Zenodo 归档拿 DOI（推荐 GitHub 集成）

1. 登录 Zenodo → Settings → GitHub → 打开该仓库的开关；
2. 回 GitHub 创建一个 Release（例如 tag `v1.0.0`），Zenodo 自动归档并分配 DOI；
3. 在该页面可拿到 DOI 徽标与 concept DOI。

完整敲除距离矩阵、20 个稳定性重复对象、CellOracle 结果不随 git 走，两种处理方式：

- **方式 A（推荐）**：这些数据可由 `code/` 重新生成，`data/README.md` 已写明输入数据来源；
  `Data and code availability` 按"代码在 GitHub、原始数据在 GEO"表述（当前 `docs/09` 即如此）；
- **方式 B**：想连派生数据一起归档，先打包再上传 Zenodo：

  ```powershell
  powershell -File code\00_setup\make_transfer_bundle.ps1 -Dest D:\transfer -Full
  Compress-Archive -Path D:\transfer\vko-project -DestinationPath D:\transfer\vko-project-full.zip
  ```

  然后把该 zip 作为 Zenodo 的一个版本上传（GitHub 单文件上限 100 MB，不要放进 git）。

## 4. 拿到 DOI 之后回填三处

1. `docs/09_methods_and_limitations_draft_EN.md` 的 `Data and code availability`：
   把 `**[待补: DOI]**` 换成真实 DOI；
2. `CITATION.cff` 增加 `doi: "10.5281/zenodo.XXXXXXX"`，并补上作者、`repository-code`、
   `preferred-citation`；
3. `README.md` 的 Citation 段可加 DOI 徽标。

## 7. 目标期刊建议（2026-09-14）

按"这是对一个已发表工具的方法学刻画 + 可复现性评估 + 可操作建议"来定位，推荐顺序：

| 顺序 | 期刊 | 理由 / 代价 |
|---|---|---|
| 1 | **PLOS Computational Biology** | 明确欢迎方法评估与可复现性研究；有真实疾病数据实例加分；开放获取，需付 APC |
| 2 | **Briefings in Bioinformatics** | 读者正是"要用这类工具的人"；适合"机制 + 三条实践建议"的写法；同样有 APC |
| 3 | **NAR Genomics and Bioinformatics** / **GigaScience** | 偏工具评估与可复现性，命中率高、竞争小，适合作为保底 |
| 4 | **Patterns**（Cell Press） | 主题最对口（scTenifoldKnk 原文即发于此），但需要按 Commentary/Perspective 改写 |
| 5 | **Genome Biology** / **Nature Methods** | 影响力最大、风险最高；Nature Methods 只适合 Correspondence/Commentary 体量 |

投稿前请核对各自当前的 scope、文章类型与 APC（期刊政策会变）。若需要，我可以联网查这几家的
现行 guide for authors 再确认一次。

## 5. 归档前自查

- [x] `code/` 与实测一致：本轮只**新增** `06_figures/54_figure_stability.R` 与
      `07_supplementary/67_stability_checkpoints_C3C4.R`，既有脚本一行未改；
- [x] 新增脚本可复现：`67_…` 归档后重跑，两个 checkpoint CSV 与归档前**逐字节一致**；
- [x] `results/figures/` 含 Figure 1–3、FigS1、**FigS2**；`results/tables/` 含 S1–**S11**；
- [x] README 的计数、图表范围、复现步骤、已知缺口已与实测一致（38 网络 / 约 20,800 次敲除）；
- [x] LICENSE：代码 MIT、`results/` 下图表 CC BY 4.0、第三方数据来源已列明；
- [ ] `CITATION.cff` 的 TBD（作者、repo、preferred-citation）——需要第 0 节的信息；
- [ ] Zenodo DOI——需要第 3 节操作；
- [ ] 若期刊要求"代码可运行性"证明，可把 `docs/12_runtime_benchmarks_CN.md` 与
      `vko-run/COMPLETION_REPORT_CN.md` 第 6 节的本机实测合并成一段补充说明。

发布前顺手修掉的三处问题（本轮已改）：

- `logs/`（本机稳定性运行日志、含绝对路径）加入 `.gitignore`，不再入库；
- `docs/08_cross_dataset_replication_CN.md` 里指向原机器绝对路径的图片链接改为相对路径
  `../results/figures/Fig_replication.png`；
- `env/python_requirements.txt` 顶部加了说明：该文件里多数条目是原机器的本地 wheel 路径，
  不能直接 `pip install -r`，版本号在 wheel 文件名里。

## 6. 本轮（2026-09-14）已完成

| 项目 | 位置 |
|---|---|
| C3/C4 检查点脚本归档 | `code/07_supplementary/67_stability_checkpoints_C3C4.R` |
| 稳定性四联图 | `code/06_figures/54_figure_stability.R` → `results/figures/FigS2_stability.{png,pdf}` |
| 稳定性逐重复表与判定表 | `results/tables/Table_S10_stability_replicates.csv`、`Table_S11_stability_summary.csv` |
| 深度匹配结果 | `outputs/results/vko_ms_depthmatched/`（27 个文件，含 `panel_depthmatched_vs_main.csv`） |
| 稿件 Methods/Results/Limitations 改写 | `docs/09_methods_and_limitations_draft_EN.md`（原稿备份在 `vko-run/backup/`） |
| 运行记录与运维脚本 | `vko-run/`（`COMPLETION_REPORT_CN.md`、`README_operations_CN.md`、日志） |
