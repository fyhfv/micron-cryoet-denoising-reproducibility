# 本地仓库已就绪：推送到 GitHub 并获取 Zenodo DOI

本目录 `reproducibility_repo` 已包含复现材料与说明页（README）。  
机器上尚未安装 `gh` 命令行，请用网页完成创建与推送（约 10 分钟）。

## A. 在 GitHub 网页创建空仓库

1. 打开 https://github.com/new  
2. Repository name 建议：`micron-cryoet-denoising-reproducibility`  
3. 选 **Public**（接 Zenodo 一般需要公开）  
4. **不要**勾选 Add README / .gitignore / license（本地已有）  
5. 点击 Create repository  
6. 复制仓库地址，例如：  
   `https://github.com/你的用户名/micron-cryoet-denoising-reproducibility.git`

## B. 把本地仓库推上去（PowerShell）

本地已完成 `git init` 与首次提交（122 个文件，约 196 MB）。你只需：

```powershell
cd "C:\Users\admin\Desktop\wy小论文\reproducibility_repo"

git remote add origin https://github.com/你的用户名/micron-cryoet-denoising-reproducibility.git
git push -u origin main
```

推送约 200 MB，需已登录 GitHub（浏览器或 Git Credential Manager）。

推送后把 `CITATION.cff` 里的 `REPLACE_ME` 改成真实用户名/仓库名，再提交一次。

## C. 连接 Zenodo（拿 DOI）

1. 打开 https://zenodo.org 并登录（可用 GitHub / ORCID）  
2. 右上角账户 → **Settings** → **GitHub**  
3. 授权后，找到本仓库，打开右侧开关（Enable）  
4. 回到 GitHub 仓库 → **Releases** → **Draft a new release**  
   - Tag：`v1.0.0`  
   - Title：`v1.0.0 — Peng_extension reproducibility`  
   - 发布说明可写：half-set recheck via `run_halves.py evaluate`；not a new algorithm  
5. 点击 **Publish release**  
6. 几分钟后回 Zenodo，打开该记录，复制 **版本 DOI**（形如 `10.5281/zenodo.xxxxxxx`）

## D. 回填正文（必做）

拿到真实 DOI 后：

1. 改本仓库 `README.md` 顶部的 Zenodo DOI 行  
2. 改稿件 *Data and code availability*：  
   - 删掉 “A permanent public archive identifier has not yet been assigned.”  
   - 把 “accompanying reproducibility files” 改成带 `https://doi.org/10.5281/zenodo.xxxxxxx` 的具体表述  

**在 DOI 出现之前不要编造链接。**

## E. 双盲审稿（若期刊要求）

可先用 Zenodo 受限访问，或录用后再公开正式 DOI；匿名阶段不要在正文写可识别个人主页的 GitHub 用户名（可用匿名仓库或仅写 DOI）。
