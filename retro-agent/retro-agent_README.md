# SimpRetro Retrosynthesis Agent

这是整个 `local_retro` 项目的根目录说明。
如果你想"启用这个 agent"，应该从当前根目录启动，而不是直接进入 `SimpRetro4Learn/` 子目录运行算法脚本。

这个 agent 的完整流程是：

1. 用户输入自然语言问题
2. LLM 解析目标分子、偏好原料和逆合成步骤数
3. agent 调用逆合成引擎
4. agent 生成自然语言解释
5. agent 输出结构图片和结果文件
6. （通过 Web 前端）可视化展示合成路线，分子结构图与 SMILES 并排显示

## 1. 你要先准备什么

在启动前，请确认你具备以下条件：

- Python 3.9
- 可以正常安装 Python 依赖
- 一个可用的 OpenAI API Key（或兼容的 API 端点，如 DeepSeek）
- 当前工作目录是项目根目录 `local_retro/`

推荐目录结构如下：

```text
local_retro/
├── README.md
├── README_EN.md
├── SKILL.md
├── requirements.txt
├── agent_cli.py
├── agent_api.py
├── agent_runtime.py
├── chem_resolution.py
├── agent_rendering.py
├── agent_config.ps1
├── run_agent.ps1
├── run_agent_api.ps1
├── smiles_to_image.py
├── static/
│   └── index.html          ← Web 前端页面
├── user_output/
│   └── agent_runs/         ← 每次运行的结果
└── SimpRetro4Learn/        ← 逆合成引擎（模板匹配 + 评分）
    ├── main.py
    ├── retro_engine.py     ← 引擎常量（数据库名、权重）
    ├── route_planner.py    ← 路线规划器（单步 + 多步）
    ├── name2smiles.py      ← 分子名称 → SMILES 解析
    ├── reaction_template.json
    ├── template_condition.json
    ├── emol_under_*.txt
    └── template/
```

## 2. 创建运行环境

推荐使用 `conda`：

```powershell
conda create -n retro_agent python=3.9
conda activate retro_agent
```

如果你不用 `conda`，也可以使用 `venv`：

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

## 3. 安装依赖

在项目根目录执行：

```powershell
pip install -r requirements.txt
```

根目录的 `requirements.txt` 会自动包含：

- root agent 层依赖
- `SimpRetro4Learn/` 中的逆合成依赖

如果安装 `rdkit` 或 `rdchiral` 失败，先解决这些底层依赖，再重新执行安装命令。这个 agent 只有在逆合成引擎依赖装齐后才能正常运行。

## 4. 配置 LLM 环境变量

你当前环境是 PowerShell，建议这样设置：

```powershell
$env:OPENAI_API_KEY="your_openai_api_key"
$env:OPENAI_MODEL="gpt-5.4-mini"
```

如果你使用兼容 OpenAI API 的代理服务，还可以额外设置：

```powershell
$env:OPENAI_BASE_URL="https://your-compatible-endpoint"
```

你也可以直接编辑 `agent_config.ps1` 填入你的配置，然后通过 `run_agent.ps1` / `run_agent_api.ps1` 一键启动。

说明：

- `OPENAI_API_KEY`：必需
- `OPENAI_MODEL`：可选，不写则使用代码中的默认模型
- `OPENAI_BASE_URL`：只有在你使用自定义网关时才需要

## 5. 启用 agent 的三种方式

你可以用三种方式启用这个 agent：

- 方式 A：命令行直接对话
- 方式 B：启动 HTTP API
- 方式 C：启动 HTTP API 后通过 Web 前端可视化使用（推荐）

### 方式 A：命令行启动

这是最简单的启用方式。

在根目录执行：

```powershell
python agent_cli.py -q "请对目标分子 CC(=O)C=C(C)C 做单步逆合成"
```

如果你希望 agent 做多步逆合成：

```powershell
python agent_cli.py -q "请对目标分子 CC(=O)C=C(C)C 做 3 步逆合成"
```

如果你希望它优先使用某些原料：

```powershell
python agent_cli.py -q "请对目标分子 CC(=O)C=C(C)C 做 2 步逆合成，优先考虑乙醇和乙酸乙酯作为起始原料"
```

命令行运行后，agent 会自动完成：

1. 识别用户语言
2. 抽取目标分子
3. 抽取逆合成步数
4. 抽取偏好原料
5. 调用逆合成规划器
6. 生成自然语言说明
7. 输出结构图片到 `user_output/agent_runs/...`

### 方式 B：启动 HTTP API

如果你准备把它接到网页、聊天界面或其他程序中，应该启动 API。

在根目录执行：

```powershell
uvicorn agent_api:app --host 0.0.0.0 --port 8000
```

启动成功后，可用接口包括：

- `GET /` — Web 前端页面
- `GET /health` — 健康检查
- `POST /agent/query` — 发送自然语言逆合成请求

你可以先检查服务是否正常：

```powershell
Invoke-RestMethod -Uri "http://127.0.0.1:8000/health"
```

然后发送自然语言请求：

```powershell
$body = @{
  message = "请对目标分子 CC(=O)C=C(C)C 做 2 步逆合成，优先使用 CCO 作为原料"
} | ConvertTo-Json

Invoke-RestMethod `
  -Uri "http://127.0.0.1:8000/agent/query" `
  -Method Post `
  -ContentType "application/json" `
  -Body $body
```

### 方式 C：Web 前端可视化（推荐）

启动 API 后，浏览器打开 `http://127.0.0.1:8000` 即可使用可视化界面。

前端功能：

- 中英文自然语言查询输入
- 步数选择器（1 步返回 3 条路线，多步返回最优路线）
- 快速示例一键填入
- 路线以**合成方向**展示（原料 → 产物），反应条件写在箭头上方，反应类型标签写在箭头下方
- 每条路线内原料分子结构图（SVG/PNG）与 SMILES 并排显示
- Score 评分、库存状态、偏好原料命中情况一目了然
- 最佳路线（排名第 1）绿色边框高亮
- 反应类型自动分类着色：Reduction / Elimination / Oxidation / Halogenation / Cycloaddition / Coupling

## 6. agent 默认行为规则

这部分很重要，决定了用户提问后 agent 的执行方式。

### 步骤数规则

- 如果用户没有明确说几步，默认按 `1` 步逆合成处理
- 如果用户要求 `1` 步逆合成，返回最推荐的 `3` 条路线
- 如果用户要求 `2` 步及以上逆合成，agent 会重复调用单步逆合成程序
- 多步模式下，最终只返回累计评分最高的 `1` 条推荐路线

### 原料偏好规则

- 如果用户给了目标原料，agent 会优先选择命中这些原料的路线
- 这不是硬过滤，而是优先排序
- 如果用户给的原料名称无法可靠解析，agent 会要求用户补充 SMILES

### 分子表达规则

- 最稳妥的输入方式是直接给 SMILES
- 如果用户给的是分子名称，agent 会尽量解析（通过 cirpy + pubchempy）
- 如果名称解析失败，agent 会要求用户明确提供结构

## 7. 推荐的用户输入方式

推荐直接这样提问：

```text
请对目标分子 CC(=O)C=C(C)C 做单步逆合成
```

或者：

```text
请对目标分子 CC(=O)C=C(C)C 做 3 步逆合成，优先使用乙醇作为原料
```

或者英文：

```text
Run a 2-step retrosynthesis for target molecule CC(=O)C=C(C)C and prefer ethanol as a starting material.
```

如果你只写"帮我逆合成阿司匹林"，agent 也许能解析名字，但这不如直接给 SMILES 稳定。

## 8. 运行后你会得到什么

每次运行结束后，结果会落盘到：

```text
user_output/agent_runs/<timestamp>_<id>/
```

通常会包含这些文件：

- `parsed_request.json`：LLM 解析出的结构化请求
- `resolved_request.json`：最终可执行的请求参数
- `planning_result.json`：逆合成规划结果
- `agent_reply.md`：给用户看的自然语言说明
- `agent_result.json`：完整输出
- `target.png` / `target.svg`：目标分子结构图（PNG + SVG）
- `route_x/` 或 `step_x/`：每条路线/每步的原料结构图（PNG + SVG）
- `final_leaf_reactants/`：多步模式下的最终原料结构图

也就是说，这个 agent 不只是"打印一段文字"，而是会同时生成说明和结构图，前端页面还可以将这些结构图嵌入合成路线图中一并展示。

## 9. 什么时候说明 agent 已经启用成功

满足下面任意一种情况，就说明 agent 已经启用了：

- `python agent_cli.py -q "..."` 能返回自然语言解释
- `uvicorn agent_api:app ...` 能正常启动
- `GET /health` 返回 `{"status": "ok"}`
- `POST /agent/query` 能返回 `reply_markdown`
- `user_output/agent_runs/` 下出现新的输出目录
- 浏览器访问 `http://127.0.0.1:8000` 能看到 Web 前端页面

## 10. 常见问题排查

### 1. 提示 `OPENAI_API_KEY is not set`

说明你还没有在当前 PowerShell 会话中设置 API Key。重新执行：

```powershell
$env:OPENAI_API_KEY="your_openai_api_key"
```

### 2. 提示缺少 `rdchiral`、`rdkit` 或其他依赖

说明底层逆合成引擎还没装好。先解决依赖安装，再运行 agent。

### 3. 用户给的是分子名称，但 agent 无法解析

这时请直接提供目标分子的 SMILES。
对于化学结构任务，SMILES 比自然语言名称更可靠。

### 4. API 能启动，但结果为空

常见原因有：

- 目标分子结构不合法
- 模板库没有覆盖到该分子
- 偏好原料过强，导致高分路线不容易命中

### 5. 为什么多步逆合成只返回一条路线

这是当前 agent 的设计规则：
多步模式下，系统会搜索多条候选，但最终只返回累计评分最高的一条最推荐路线。

## 11. 你真正需要运行的入口

如果你的目标是"启用整个 agent"，你真正应该用的是根目录入口：

- `agent_cli.py` — 命令行
- `agent_api.py` — HTTP API + Web 前端

而不是直接运行子目录里的：

- `SimpRetro4Learn/main.py`

后者更偏底层算法测试；前者才是完整的 LLM agent 工作流入口。

## 12. 当前限制

- 底层模板引擎本质上仍是单步逆合成
- 多步能力来自上层规划搜索（迭代调用单步引擎）
- 输出是启发式建议，不是实验验证方案
- 如果目标结构无法可靠解析，仍然需要用户补充 SMILES

## 13. 最短启用步骤

如果你只想最快跑起来，按下面 5 步做：

```powershell
conda create -n retro_agent python=3.9
conda activate retro_agent
pip install -r requirements.txt
$env:OPENAI_API_KEY="your_openai_api_key"
python agent_cli.py -q "请对目标分子 CC(=O)C=C(C)C 做单步逆合成"
```

要使用 Web 可视化界面，再加一步：

```powershell
python -m uvicorn agent_api:app --host 127.0.0.1 --port 8000
# 浏览器打开 http://127.0.0.1:8000
```

做到这里，这个 agent 就已经启用了。
