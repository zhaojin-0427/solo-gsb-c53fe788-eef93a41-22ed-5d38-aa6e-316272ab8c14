# 卫星过境观测 API（Satellite Pass Observation API）

给定一个历元时刻的卫星 ECI 状态与地面观测站坐标，预测指定 UTC 时间区间内
卫星相对观测站的可见窗口（升/落时刻、峰值仰角及其时刻、方位角等）。

- **语言/框架**：Python 3.11+、FastAPI、NumPy、Pydantic v2
- **动力学**：纯两体开普勒传播（WGS84 引力常数），牛顿迭代解开普勒方程
  （带回退二分法），机器精度守恒
- **坐标系统**：固定 WGS84 椭球与固定地球自转率
  （ω = 7.2921151467×10⁻⁵ rad/s），IAU 1982 GMST 模型做 ECI↔ECEF 旋转；
  观测站按 WGS84 大地坐标转换为 ECEF，再到站心地平坐标（东-北-天，ENU）
- **求交策略**：两级自适应扫描 + 二分精确求仰角穿越 + 黄金分割精化峰值，
  搜索区间按**半开区间 [start, end)** 处理；双端低值区间须经 Lipschitz
  变化率界证明不可见才计为已覆盖，否则递归细分寻找隐藏穿越，细分到
  1 µs 仍不可判定的擦边区间如实登记 `uncovered_intervals`
  （见下文“算法说明”）

---

## 1. 一键启动（Docker Compose）

```bash
docker compose up --build
```

服务监听 `http://localhost:8000`：

| 用途 | 地址 |
| --- | --- |
| 健康检查 | `GET http://localhost:8000/health` |
| 过境预测（单站） | `POST http://localhost:8000/passes` |
| 接力计划（多站） | `POST http://localhost:8000/relay-plan` |
| 交互式文档（Swagger UI） | `http://localhost:8000/docs` |
| OpenAPI 模式 | `http://localhost:8000/openapi.json` |

停止服务：`docker compose down`。

### 不使用 Docker（本地运行）

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

运行测试：

```bash
pip install -r requirements.txt pytest httpx
pytest -q
```

---

## 2. 请求格式（POST /passes）

Content-Type: `application/json`。所有时间必须是**带时区的 UTC**
（ISO-8601，推荐以 `Z` 结尾）。

| 字段 | 类型 | 单位 / 范围 | 说明 |
| --- | --- | --- | --- |
| `epoch` | string | UTC 时间 | 卫星状态向量的历元 |
| `start` / `end` | string | UTC 时间 | 搜索区间，`end > start`，跨度 **≤ 7 天** |
| `min_elevation_deg` | number | 度，(-90, 90) | 最低仰角（含该仰角即为可见） |
| `satellite.position_km` | number[3] | **km**，有限数 | 历元 ECI 位置 [x, y, z] |
| `satellite.velocity_km_s` | number[3] | **km/s**，有限数 | 历元 ECI 速度 [vx, vy, vz] |
| `station.latitude_deg` | number | 度，[-90, 90] | WGS84 大地纬度 |
| `station.longitude_deg` | number | 度，[-180, 180] | WGS84 经度 |
| `station.altitude_km` | number | **km**，[-12, 100000] | WGS84 椭球高 |

### 422 校验规则

以下情况返回 **HTTP 422**（FastAPI/Pydantic 标准错误体）：

- 时间顺序非法（`end <= start`）或跨度超过 7 天；时间不带时区
- 经纬高、最低仰角超出物理范围；向量长度不是 3 或含 NaN/Infinity
- 轨道退化：位置/速度为零、r 与 v 共线（角动量为零）
- 轨道非束缚：抛物线/双曲线（比机械能 ≥ 0）或离心率 ≥ 1
- 近地点位于地表以下（两体传播无物理意义）
- 未知字段（请求体禁止额外字段）

> 单位提示：位置必须用 **km**、速度必须用 **km/s**。若误把位置按米
> （放大 1000 倍）提交，近地点将落到地球内部，服务会以 422
> `orbit_impact_earth` 拒绝；速度按 m/s 提交则会因双曲线轨道被拒绝。

---

## 3. 实际请求示例

以下是一次真实请求（420 km 圆轨道、倾角 51.6° 的类 ISS 卫星，
北京附近观测站，最低仰角 10°，搜索 6 小时；完整文件见
[`examples/request.json`](examples/request.json)）：

```bash
curl -X POST http://localhost:8000/passes \
  -H 'Content-Type: application/json' \
  -d @examples/request.json
```

请求体（即 `examples/request.json`）：

```json
{
  "epoch": "2025-06-01T00:00:00Z",
  "start": "2025-06-01T00:00:00Z",
  "end": "2025-06-01T06:00:00.000Z",
  "min_elevation_deg": 10.0,
  "satellite": {
    "position_km": [4270.307144917519, -1323.5276198332833, 5121.271150639637],
    "velocity_km_s": [3.950957000514434, 6.347260232247907, -1.6540865405752112]
  },
  "station": {
    "latitude_deg": 39.9,
    "longitude_deg": 116.4,
    "altitude_km": 0.05
  }
}
```

### 实际响应（HTTP 200）

与 [`examples/response.json`](examples/response.json) 内容一致：

```json
{
  "status": "complete",
  "count": 1,
  "windows": [
    {
      "start": "2025-06-01T00:01:40.294555Z",
      "end": "2025-06-01T00:08:29.127935Z",
      "peak_time": "2025-06-01T00:05:05.712339Z",
      "peak_elevation_deg": 85.67588427634527,
      "peak_azimuth_deg": 218.30123915031336,
      "rise_azimuth_deg": 306.5090307458435,
      "set_azimuth_deg": 130.01647942552455,
      "max_range_km": 1521.862803277411,
      "duration_s": 408.8333794882783,
      "start_open": false,
      "end_open": false,
      "truncated": false
    }
  ],
  "uncovered_intervals": [],
  "search_start": "2025-06-01T00:00:00Z",
  "search_end": "2025-06-01T06:00:00Z",
  "min_elevation_deg": 10.0,
  "diagnostics": {
    "orbit": {
      "semi_major_axis_km": 6798.136999999997,
      "eccentricity": 3.6354953649115305e-16,
      "period_s": 5578.222707272108,
      "perigee_radius_km": 6798.136999999995,
      "apogee_radius_km": 6798.137,
      "perigee_altitude_km": 419.99999999999545,
      "apogee_altitude_km": 420.0
    },
    "coarse_samples": 297,
    "dense_samples": 1898,
    "search_span_s": 21600.0,
    "uncovered_seconds": 0,
    "coverage_fraction": 1.0,
    "notes": "Elevation crossings refined by bisection; peaks by golden-section search. Half-open search interval [start, end). Below-threshold intervals count as covered only when a Lipschitz rate bound proves them invisible; unprovable brackets are bisected down to the 1 us reporting resolution to expose hidden crossings, and whatever stays unresolvable (grazing contacts) is reported in uncovered_intervals instead of being silently skipped."
  }
}
```

字段含义：

- `windows[]` 按 **`start` 升序**排列；每个窗口含起止时刻、峰值时刻、
  峰值仰角（度）、峰值/升/落方位角（度，正北顺时针）、最大斜距（km）、
  持续时间（秒）。
- `start_open` / `end_open`：窗口在搜索区间起点已可见 / 到区间终点仍可见
  （半开区间被边界截断）。
- `truncated`：窗口因数值失败被提前截断（此时 `status` 为 `partial`）。
- `uncovered_intervals`：未能覆盖的时间区间及原因——传播/求根数值失败，
  或**无法在 1 µs 报告精度内判定**的擦边接触（峰值与门限之差小于数值
  分辨率，既不能证明不可见也无法解析出穿越）。不为空时 `status` 为
  `partial`（有结果但不完整）或 `failed`（没有任何可信窗口），
  `coverage_fraction` 相应小于 1。
- `diagnostics`：轨道根数、粗/细扫描采样数、未覆盖秒数与覆盖率。

### 422 响应示例（搜索区间 8 天 > 7 天上限）

将 `examples/request.json` 的 `end` 改为 `"2025-06-09T00:00:00Z"` 后
（完整响应见 [`examples/error_422.json`](examples/error_422.json)）：

```json
{
  "detail": [
    {
      "type": "value_error",
      "loc": ["body"],
      "msg": "Value error, search span (end - start) must not exceed 7 days",
      "input": { "epoch": "2025-06-01T00:00:00Z", "...": "..." },
      "ctx": { "error": {} }
    }
  ]
}
```

轨道退化/非束缚错误示例：

```json
{
  "detail": [
    {
      "type": "orbit_degenerate_or_invalid",
      "loc": ["body", "satellite"],
      "msg": "orbit is not bound (parabolic/hyperbolic): only elliptic orbits are supported"
    }
  ]
}
```

---

## 4. 多站接力计划（POST /relay-plan）

给定与 `/passes` 相同的历元、ECI 卫星状态与 UTC 搜索区间，以及一个
**非空的测站列表**（每站带**唯一 ID**、坐标和本站的最低仰角），服务对每站
**复用同一套两体传播与 `find_passes`** 求可见性，然后在所有站“已证可见”的
窗口上规划搜索期内的**最长连续接力时段**：由哪些站各负责哪一段、在什么时刻
交接。

### 4.1 接力规则（只使用已证可见）

- 每个站的可见性**逐站**调用现有的两体传播 + `find_passes` 得到。站级
  `status`、**半开**窗口 `[start, end)`、`truncated` 截断标志与
  `uncovered_intervals` 一律保留，窗口本身即“已证可见时段”。
- **交接只能在两个站窗口有正时长重叠时发生**：
  `max(start_A, start_B) < min(end_A, end_B)`。仅在某一瞬时相接的两个窗口
  不能交接（半开区间下该点不属于前一站）。切换时刻取前一站窗口的 `end`，
  该时刻严格落在后一站窗口内部，因此交接无缝、双方在切换前后均已证可见。
- 目标是搜索期内的**最长连续接力时段**；等长时先选**切换次数少**的方案，
  再按**开始时刻更早**、负责站 **ID 元组字典序最小**稳定裁决。
- 某站自己无法覆盖、但**别的站在该时段已证可见**时，接力照常进行（该时段
  不算未知）。
- **无人可证可见、但至少有一站在该时段无法判定**（出现在其
  `uncovered_intervals` 中）时，该段列入 `unknown_intervals`（附各站原因），
  `status` 如实降级——**既不报“确定不可见”，也不报“完整覆盖”**。所有站都
  已用变化率界证明低于门限的时段才是确定不可见，不列入未知区间。
- 没有任何可见窗口时返回**空计划**（`plan: null`）：同时存在未判定时段则
  `status` 为 `failed`，全部确定不可见则 `status` 为 `complete`。

### 4.2 请求格式

Content-Type: `application/json`，所有时间为带时区 UTC。

| 字段 | 类型 | 单位 / 范围 | 说明 |
| --- | --- | --- | --- |
| `epoch` | string | UTC 时间 | 卫星状态向量的历元（同 `/passes`） |
| `start` / `end` | string | UTC 时间 | 搜索区间，`end > start`，跨度 **≤ 7 天** |
| `satellite.position_km` | number[3] | km，有限数 | 历元 ECI 位置 |
| `satellite.velocity_km_s` | number[3] | km/s，有限数 | 历元 ECI 速度 |
| `stations` | array | **非空** | 候选接力站列表 |
| `stations[].id` | string | 非空、**唯一** | 站 ID，用于计划与交接标注 |
| `stations[].latitude_deg` | number | 度，[-90, 90] | WGS84 大地纬度 |
| `stations[].longitude_deg` | number | 度，[-180, 180] | WGS84 经度 |
| `stations[].altitude_km` | number | km，[-12, 100000] | WGS84 椭球高 |
| `stations[].min_elevation_deg` | number | 度，(-90, 90) | **该站**的最低仰角（含该仰角即可见） |

**422 规则**：与 `/passes` 相同的时间/坐标/向量/轨道校验（退化、非束缚、
近地点入地均 422），另外 `stations` 为空、站 `id` 为空/纯空白、或站 `id`
重复也返回 422。

### 4.3 实际请求示例

类 ISS 的 420 km 圆轨道，北京（`BEIJING`，116.4°E）与东京湾
（`TOKYO-BAY`，139.7°E）两站，最低仰角均 10°，搜索区间取两段过境的并集
（完整文件见 [`examples/relay_request.json`](examples/relay_request.json)）：

```bash
curl -X POST http://localhost:8000/relay-plan \
  -H 'Content-Type: application/json' \
  -d @examples/relay_request.json
```

请求体（即 `examples/relay_request.json`）：

```json
{
  "epoch": "2025-06-01T00:00:00Z",
  "start": "2025-06-01T00:01:40.294555Z",
  "end": "2025-06-01T00:09:41.322998Z",
  "satellite": {
    "position_km": [4270.307144917519, -1323.5276198332833, 5121.271150639637],
    "velocity_km_s": [3.950957000514434, 6.347260232247907, -1.6540865405752112]
  },
  "stations": [
    {"id": "BEIJING",   "latitude_deg": 39.9, "longitude_deg": 116.4,
     "altitude_km": 0.05, "min_elevation_deg": 10.0},
    {"id": "TOKYO-BAY", "latitude_deg": 35.5, "longitude_deg": 139.7,
     "altitude_km": 0.02, "min_elevation_deg": 10.0}
  ]
}
```

### 4.4 实际响应（HTTP 200）

与 [`examples/relay_response.json`](examples/relay_response.json) 一致：

```json
{
  "status": "complete",
  "search_start": "2025-06-01T00:01:40.294555Z",
  "search_end": "2025-06-01T00:09:41.322998Z",
  "plan": {
    "start": "2025-06-01T00:01:40.294555Z",
    "end": "2025-06-01T00:09:41.322998Z",
    "duration_s": 481.028443,
    "switches": 1,
    "assignments": [
      {
        "station_id": "BEIJING",
        "start": "2025-06-01T00:01:40.294555Z",
        "end": "2025-06-01T00:08:29.127935Z",
        "duration_s": 408.83338,
        "truncated": false
      },
      {
        "station_id": "TOKYO-BAY",
        "start": "2025-06-01T00:08:29.127935Z",
        "end": "2025-06-01T00:09:41.322998Z",
        "duration_s": 72.195063,
        "truncated": false
      }
    ],
    "handoffs": [
      {
        "time": "2025-06-01T00:08:29.127935Z",
        "from_station_id": "BEIJING",
        "to_station_id": "TOKYO-BAY"
      }
    ],
    "stations": ["BEIJING", "TOKYO-BAY"]
  },
  "stations": [
    {
      "id": "BEIJING",
      "status": "complete",
      "truncated": false,
      "windows": [
        {"start": "2025-06-01T00:01:40.294555Z",
         "end": "2025-06-01T00:08:29.127935Z", "truncated": false}
      ],
      "uncovered_intervals": []
    },
    {
      "id": "TOKYO-BAY",
      "status": "complete",
      "truncated": false,
      "windows": [
        {"start": "2025-06-01T00:07:28.047544Z",
         "end": "2025-06-01T00:09:41.322998Z", "truncated": false}
      ],
      "uncovered_intervals": []
    }
  ],
  "unknown_intervals": []
}
```

**如何读取切换时刻**：两站窗口在 `00:07:28`–`00:08:29` 正时长重叠
（约 61 秒），计划在 `plan.handoffs[0].time = 00:08:29.127935Z` 把责任从
`BEIJING` 交给 `TOKYO-BAY`——这正是北京窗口的结束时刻，且严格落在东京湾
窗口内部。`plan.assignments[]` 是无缝拼接的负责区间（前一段 `end` 等于
后一段 `start`），每段给出负责站与该段时长；`plan.switches` 是切换次数，
`plan.stations` 是负责站序列。本例接力覆盖整个搜索区间且无未知时段，
故顶层 `status` 为 `complete`；搜索区间比接力段长时为 `partial`。

其余字段：

- `stations[]`：按请求顺序返回**每站**的 `status`（来自 `find_passes`）、
  半开已证可见窗口（含 `truncated`）与该站自己的 `uncovered_intervals`。
- `unknown_intervals[]`：无人已证可见、但至少有一站无法判定的时段及各站
  原因（`"站ID: 原因"`）。所有站均证明低于门限的确定不可见时段不在其中。
- 没有任何可见窗口时 `plan` 为 `null`：有未知时段则 `status="failed"`，
  全部确定不可见则 `status="complete"`。

### 4.5 接力 422 响应示例（站 ID 重复）

把第二个站的 `id` 也改成 `"BEIJING"`（完整响应见
[`examples/relay_error_422.json`](examples/relay_error_422.json)）：

```json
{
  "detail": [
    {
      "type": "value_error",
      "loc": ["body"],
      "msg": "Value error, station 'id' values must be unique; duplicates: BEIJING"
    }
  ]
}
```

---

## 5. 算法说明

### 5.1 两体传播

由历元状态计算比机械能 ε = v²/2 − μ/r 与离心率矢量，要求 ε < 0、
e < 1。任意时刻 t = t₀ + Δt 用开普勒方程 M = E − e·sin E
（牛顿迭代，容差 1e-10，失败时二分回退），再在由 h、e 矢量构造的
近焦点坐标系中还原 r、v。圆轨道（e≈0）使用 r 方向作为面内基准，
无近点奇点。自测：一个周期后状态误差为 0，3 天内能量/角动量
相对漂移 < 1e-15。

### 5.2 坐标与时角

- GMST 使用 IAU 1982 模型（J2000 锚点 18h41m50.543s）。距 J2000 的
  日数直接由秒差计算——若先合成完整儒略日（≈2.46×10⁶）再相减，双精度
  会把时刻量化到约 40 µs，经测站杠杆臂（R站/斜距，近天顶时可达 12）
  放大为约 1e-6° 的仰角阶梯噪声，足以淹没或伪造擦边穿越；
- ECI→ECEF：`r_ECEF = Rz(−GMST)·r_ECI`，速度另减 ω×r 项；
- ENU 基由测站纬经度固定，仰角 `asin(ρ_up/|ρ|)`、
  方位角 `atan2(ρ_E, ρ_N)`；并给出仰角/方位角解析变化率供自适应用。

### 5.3 两级自适应扫描（不漏窗口）

1. **粗扫（几何门限）**：只算卫星—测站的**地心角** α。可见要求
   cos α ≥ R站/r星，因此用远地点构造最松门限 ψ_flag（地平角 + 8° 缓冲）。
   由近地点速度给出严格的地心角变化率上界
   Kα = v_p/r_p + ω_E，粗扫步长固定为 5°/Kα（LEO 约 40–75 s），
   **任何可能可见的几何都不可能被跨过**。标记点 ±一个粗步长合并为候选段。
2. **细扫（仰角门限）**：在候选段内按仰角距阈值的余量和全局 Lipschitz 界
   K1 = (v_p + ω_E·r_a)/(r_p − R站) 自适应取步（下限 0.02 s、上限 20 s）。
   **双端低于阈值的区间只有满足 K1·Δt ≤ 0.5·min(余量) 才算已覆盖**
   （变化率界证明曲线不可能触及阈值）；步长下限本身**不构成**覆盖依据。
3. **不可证明区间的递归细分**：不满足证明的双端低值区间被递归二分
   ——中点若高于阈值即暴露隐藏穿越（交回二分精确化），子区间若被证明
   则覆盖，直到 1 µs（UTC 时间戳的报告精度）。在此精度下仍不可证明的
   残段（峰值与门限之差小于数值分辨率的擦边接触）**登记进
   `uncovered_intervals` 并注明原因**，绝不静默跳过；细分求值设有全局
   预算（20 万次，正常 7 天 LEO 任务约用 3 万次），预算耗尽时剩余未证明
   区间同样如实登记。
4. **精确求交**：符号变化的相邻采样对做 45 次二分（时间精度 1e-8 s，
   端点仰角误差 < 1e-7°）。
5. **峰值**：以最佳细扫采样及其相邻采样为支架做 32 次黄金分割，
   且结果永不低于实际观测到的最大采样值。

保证的边界：只有峰值与阈值之差小于数值分辨率（仰角噪声 <~1e-7°，
见 4.2）、持续时间短到无法在 1 µs 精度内表示的擦边接触才**可能**无法
判定——此时区间与原因会出现在 `uncovered_intervals` 中，`status` 与
`coverage_fraction` 如实降级，而不会以步长下限为由宣称 `complete`；
其余窗口不漏。与独立的 0.25 s 均匀网格真值对比（1 天类 ISS 轨道）：
窗口数完全一致，升/落时刻误差 < 0.05 s，峰值仰角误差 < 0.001°。

### 5.4 半开区间与失败报告

- 搜索域为 **[start, end)**：恰在 `end` 发生的降落穿越按区间内处理
  （窗口闭合于该时刻）；恰在 `start` 已可见则标记 `start_open=true`，
  到 `end` 仍可见标记 `end_open=true`。穿越时刻按 UTC 微秒精度报告；
  距 `start` 不足 1 µs 的上升穿越与恰在边界上不可区分，按
  `start_open=true` 处理。
- 零时长切触（峰值恰好等于阈值）不作为窗口返回。
- 任何传播/求值/二分失败：受影响时段进入 `uncovered_intervals`
  （含起止时刻与原因），打开中的窗口标记 `truncated`，
  `status` 降级为 `partial`/`failed`，绝不静默丢弃。
- 双端低值且无法用变化率界证明不可见、细分到 1 µs 仍无法判定的
  擦边区间同样进入 `uncovered_intervals`（相邻残段合并为一条），
  `coverage_fraction = 1 − 未覆盖秒数/搜索跨度` 随之如实下降。

### 5.5 多站接力规划（relay.py）

- 对每个站独立调用 `find_passes`，其窗口即该站的“已证可见”半开区间，
  `truncated` 与 `uncovered_intervals` 原样保留，**不做任何乐观外推**。
- 以窗口为节点、**异站且正时长重叠**为边，用并查集求连通分量；每个分量的
  `[最早 start, 最晚 end)` 即一条可连续接力的时段（贪心“最远可达”论证保证
  存在一条正重叠交接链贯穿全段）。
- 分量内按 `end` 升序做 DAG 最短路求**最少切换次数**的窗口链：边 j→i 要求
  异站、`w_i.start < w_j.end < w_i.end`（i 正时长包住 j 的结束点）且
  `w_i.end > w_j.end`（i 把接力延伸得更远）。切换点取前一窗口的 `end`，
  它严格落在后一窗口内部，故负责区间无缝拼接。
- 等长方案按 `(−时长, 切换次数, 开始时刻, 负责站 ID 元组)` 取最小稳定裁决。
- 未知时段 = 搜索区间减去所有站可见窗口并集，再与所有站
  `uncovered_intervals` 并集**取交**：只有“无人证明可见且有人无法判定”的
  部分才列入 `unknown_intervals`；所有站都已证明低于门限的部分是确定
  不可见，绝不报未知，也绝不报完整覆盖。

### 5.6 性能

典型开销（单机纯 Python）：7 天类 ISS 任务约 8300 次粗扫 + 6.7 万次
细扫求值，耗时约 7 秒；6 小时任务约 0.1–0.5 秒。高轨（GEO、闪电轨道）
更省或相当。病态擦边场景（峰值距阈值 <1e-6°）会触发递归细分，
单次请求额外求值以 20 万次为上限。`/relay-plan` 的传播/求值开销随站数
线性增长（每站一次完整的 `/passes` 级计算），接力图本身为 O(W²)。

---

## 6. 项目结构

```
app/
  main.py        # FastAPI 路由、422 映射、/passes 与 /relay-plan 响应组装
  schemas.py     # Pydantic 请求/响应模型与范围校验（含站 ID 唯一性）
  mechanics.py   # WGS84 常数、两体开普勒传播
  geometry.py    # 大地坐标、GMST、ECI/ECEF、ENU、仰角/方位角
  passes.py      # 两级自适应扫描、二分穿越、黄金分割峰值
  relay.py       # 已证可见窗口的连通分量、最少切换 DAG、未知区间
  time_utils.py  # UTC/儒略日/GMST
examples/        # /passes 与 /relay-plan 的真实请求、成功与 422 响应
tests/           # pytest：52 项（网格真值交叉验证、擦边回归、接力交接/裁决/未知语义、422 全分支）
Dockerfile, docker-compose.yml, requirements.txt
```
