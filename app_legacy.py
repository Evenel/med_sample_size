import math
import streamlit as st
import ui_continuous
import ui_proportion
import ui_survival
from calc_proportion import ProportionCalculator, ProportionInput
from calc_continuous import ContinuousCalculator
from calc_survival import SurvivalCalculator

# 1. 页面基础设置
st.set_page_config(
    page_title="临床试验样本量计算",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 2. 注入修复后的 CSS
def inject_morandi_ui():
    st.markdown(
        """
        <style>
        /* 全局背景 */
        body, .stApp {
            background-color: #F4F2ED;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            color: #3C3C3C;
        }
        
        /* 修复侧边栏不见的问题：去掉了 header 的隐藏 */
        #MainMenu, footer {visibility: hidden;}
        
        /* 居中大标题 */
        .main-title {
            font-size: 32px;
            font-weight: 700;
            color: #3C3C3C;
            margin-bottom: 40px;
            margin-top: -20px;
            letter-spacing: 1px;
            text-align: center;
        }
        
        /* 紧密相连的标题线（去掉横线，用左侧边框区分） */
        .section-title {
            font-size: 16px;
            font-weight: 600;
            color: #4A4A4A;
            margin-bottom: 16px;
            margin-top: 32px;
            border-left: 3px solid #788A9A;
            padding-left: 8px;
            line-height: 1;
        }
        .section-title:first-child { margin-top: 0; }
        
        /* 结果展示大卡片 */
        .result-kpi-card {
            background-color: #788A9A;
            border-radius: 12px;
            padding: 28px 24px;
            color: #FFFFFF;
            box-shadow: 0 8px 16px rgba(120, 138, 154, 0.2);
            margin-bottom: 24px;
        }
        
        /* 调整后的分组人数视觉平衡 */
        .group-split {
            display: flex;
            justify-content: space-between;
            gap: 16px;
            margin-top: 20px;
            margin-bottom: 24px;
        }
        .group-box {
            background: rgba(255, 255, 255, 0.15);
            border-radius: 8px;
            padding: 20px 16px;
            flex: 1;
            text-align: center;
            border: 1px solid rgba(255, 255, 255, 0.2);
        }
        .group-num { font-size: 36px; font-weight: 700; line-height: 1.1; margin-bottom: 8px; }
        .group-name { font-size: 16px; font-weight: 600; opacity: 1; letter-spacing: 0.5px; }
        
        /* 总人数汇总栏文字 */
        .total-box {
            text-align: center;
            font-size: 18px; 
            font-weight: 600;
            border-top: 1px dashed rgba(255, 255, 255, 0.4);
            padding-top: 20px;
            letter-spacing: 0.5px;
        }

        /* 专业话术框，支持段落缩进或空行 */
        .statement-box {
            background-color: #EBE7E0;
            border-radius: 8px;
            padding: 24px;
            font-size: 14px;
            line-height: 1.7;
            color: #3C3C3C;
            border: 1px solid #DFDBD3;
        }
        .statement-box p { margin-bottom: 12px; text-indent: 2em; }
        .statement-box p:last-child { margin-bottom: 0; }
        
        /* 控件样式微调 */
        div[data-testid="stNumberInput"] label p, div[data-testid="stSelectbox"] label p {
            color: #5A5A5A;
            font-size: 13px;
            font-weight: 500;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

# 3. 侧边栏
def sidebar_layout():
    with st.sidebar:
        st.markdown("<h3 style='color: #4A4A4A; font-size: 17px; margin-bottom: 20px;'>试验设计导航</h3>", unsafe_allow_html=True)
        variable_type = st.selectbox("主要终点变量", ["分类变量 (率)", "连续变量 (均值)", "生存时间 (HR)"], index=0)
        simon_design_type = None

        if variable_type == "分类变量 (率)":
            group_type = st.selectbox(
                "试验分组设计",
                ["独立两组 (Parallel)", "单臂设计 (Single Arm)", "配对两组 (Paired)", "三组设计 (Three Arm)", "Simon 二阶段"],
                index=0,
            )
            if group_type == "Simon 二阶段":
                test_type = "优效性检验"
                simon_design_type = st.selectbox(
                    "Simon 设计类型",
                    ["optimal", "minimax"],
                    index=0,
                    help="optimal: 最小化 EN(p0)；minimax: 最小化最大样本量 n。",
                )
                st.text_input("假设检验目标", value="Simon 二阶段设计", disabled=True)
            elif group_type == "三组设计 (Three Arm)":
                test_type = st.selectbox(
                    "假设检验目标",
                    [
                        "差异性检验（率齐性）",
                        "有序趋势检验（Cochran-Armitage）",
                    ],
                    index=0,
                    help="齐性：任意组间率是否不同；趋势：三组有自然顺序（如剂量）时检验线性趋势。",
                )
            else:
                if group_type in ("单臂设计 (Single Arm)", "独立两组 (Parallel)"):
                    test_type = st.selectbox(
                        "假设检验目标",
                        ["优效性检验", "非劣效性检验", "差异性检验", "等效性检验", "精度分析"],
                        index=0,
                    )
                else:
                    test_type = st.selectbox(
                        "假设检验目标",
                        ["优效性检验", "非劣效性检验", "差异性检验", "等效性检验"],
                        index=0,
                    )
        elif variable_type == "连续变量 (均值)":
            group_type = st.selectbox(
                "试验分组设计",
                ["独立两组 (Parallel)", "单臂设计 (Single Arm)", "配对两组 (Paired)", "三组设计 (Three Arm)"],
                index=0,
            )
            if group_type == "三组设计 (Three Arm)":
                test_type = st.selectbox("假设检验目标", ["差异性检验"], index=0, help="临床上三组均值常用方差分析总体差异检验。")
            else:
                test_type = st.selectbox("假设检验目标", ["优效性检验", "非劣效性检验", "差异性检验", "等效性检验"], index=0)
        else:
            group_type = st.selectbox("试验分组设计", ["独立两组 (Parallel)", "单臂设计 (Single Arm)"], index=0)
            test_type = st.selectbox("假设检验目标", ["优效性检验", "非劣效性检验"], index=0)

    return variable_type, group_type, test_type, simon_design_type

# 4. 主工作区
def main_layout(variable_type, group_type, test_type, simon_design_type):
    inject_morandi_ui()
    st.markdown('<div class="main-title">临床试验样本量计算</div>', unsafe_allow_html=True)

    supported_combos = {
        ("分类变量 (率)", "独立两组 (Parallel)"),
        ("分类变量 (率)", "单臂设计 (Single Arm)"),
        ("分类变量 (率)", "Simon 二阶段"),
        ("分类变量 (率)", "三组设计 (Three Arm)"),
        ("连续变量 (均值)", "独立两组 (Parallel)"),
        ("连续变量 (均值)", "单臂设计 (Single Arm)"),
        ("连续变量 (均值)", "配对两组 (Paired)"),
        ("分类变量 (率)", "配对两组 (Paired)"),
        ("生存时间 (HR)", "独立两组 (Parallel)"),
        ("生存时间 (HR)", "单臂设计 (Single Arm)"),
    }
    if (variable_type, group_type) not in supported_combos:
        st.warning(f"当前组合【{variable_type} + {group_type}】已纳入菜单，但算法模块仍在开发中。")
        return

    col_input, col_result = st.columns([5, 5], gap="large")

    # ------------------ 左侧动态表单区 ------------------
    # 默认值（便于后端调度）
    single_arm_method = "score"  # score / exact（仅单臂率功效分析使用）
    half_width_percent = 10.0
    precision_method = "正态近似法 (Wald)"
    two_arm_diff_method = None  # None=auto, pooled_z, fisher_exact, arcsine
    two_arm_margin_method = None  # None=auto, fm_score, wald_unpooled, unconditional_chan
    paired_rate_method = None  # None=auto, nam_score, exact_mcnemar, wald_blackwelder, exact_conditional_paired
    three_arm_diff_method = None  # None=auto, pearson, monte_carlo（仅率齐性）

    with col_input:
        # 彻底删除了错误的 <div class="morandi-card"> HTML 标签，解决空框 Bug
        st.markdown('<div class="section-title">临床预期</div>', unsafe_allow_html=True)
        
        if variable_type == "分类变量 (率)" and group_type == "Simon 二阶段":
            c_p1, c_p2 = st.columns(2)
            with c_p1:
                p_treatment_percent = st.number_input("目标有效率 P1 (%)", value=40.0, step=0.1)
            with c_p2:
                p_control_percent = st.number_input("无效界值 P0 (%)", value=20.0, step=0.1)
            margin_percent = 0.0
        elif variable_type == "连续变量 (均值)":
            _c = ui_continuous.render_clinical_inputs(group_type, test_type)
            mean_t = _c.get("mean_t")
            mean_c = _c.get("mean_c")
            sd_t = _c.get("sd_t")
            sd_c = _c.get("sd_c")
            mean_ref = _c.get("mean_ref")
            sd_single = _c.get("sd_single")
            mean_diff = _c.get("mean_diff")
            sd_diff = _c.get("sd_diff")
            margin_percent = _c["margin_percent"]
        elif variable_type == "生存时间 (HR)":
            _s = ui_survival.render_clinical_inputs(group_type, test_type, legacy_simple=True)
            median_t = _s["median_t"]
            event_rate_percent = _s["event_rate_percent"]
            margin_percent = _s["margin_percent"]
            if group_type == "单臂设计 (Single Arm)":
                median_ref = _s["median_ref"]
            else:
                median_c = _s["median_c"]
        elif variable_type == "分类变量 (率)" and group_type == "配对两组 (Paired)":
            _p = ui_proportion.render_clinical_inputs(group_type, test_type, simon_design_type)
            p10_pct = _p["p10_pct"]
            p01_pct = _p["p01_pct"]
            margin_percent = _p["margin_percent"]
            direction = _p["direction"]
        elif variable_type == "分类变量 (率)" and group_type == "三组设计 (Three Arm)":
            c3a, c3b, c3c = st.columns(3)
            with c3a:
                p1_pct = st.number_input("第 1 组有效率 π₁ (%)", value=50.0, step=0.1)
            with c3b:
                p2_pct = st.number_input("第 2 组有效率 π₂ (%)", value=40.0, step=0.1)
            with c3c:
                p3_pct = st.number_input("第 3 组有效率 π₃ (%)", value=35.0, step=0.1)
            margin_percent = 0.0
            direction = "higher"
            if test_type == "差异性检验（率齐性）":
                st.caption(
                    "三组独立二分类：率齐性检验 H₀：π₁=π₂=π₃（任意组间存在差异即可能拒绝）。等额分配；"
                    "常规用 Pearson χ² 非中心；若期望频数不足或率极偏，可用蒙特卡洛。"
                )
                with st.expander("三组率齐性检验（算法说明）", expanded=False):
                    st.markdown(
                        """
| 场景 | 推荐 | 说明 |
|------|------|------|
| 期望频数较充分 | **Pearson χ² 齐性**（非中心 χ² 样本量） | 与 2×2 表 Pearson χ² 同一渐近族 |
| 小样本 / 极端率 / 需模拟复核 | **蒙特卡洛**（Pearson 统计量 + 渐近 p） | 模拟 3×2 表估计功效，注册试验建议统计师复核 |

* **自动**：三组率均在 [20%,80%] 且按 Pearson 解得的每组最小期望频数 ≥5 → Pearson；否则 → 蒙特卡洛。*
"""
                    )
            else:
                st.caption(
                    "三组有序（如安慰剂→低剂量→高剂量）：**Cochran–Armitage** 线性趋势检验（双侧）。"
                    "分数须严格递增且与临床顺序一致；样本量由蒙特卡洛估计（渐近 Z 检验 p 值）。"
                )
                cs1, cs2, cs3 = st.columns(3)
                with cs1:
                    ca_t1 = st.number_input("第 1 组趋势分数", value=0.0, step=0.5)
                with cs2:
                    ca_t2 = st.number_input("第 2 组趋势分数", value=1.0, step=0.5)
                with cs3:
                    ca_t3 = st.number_input("第 3 组趋势分数", value=2.0, step=0.5)
                with st.expander("Cochran–Armitage 趋势检验说明", expanded=False):
                    st.markdown(
                        """
- **适用**：组别有**自然顺序**（剂量、暴露等级等），且主要关心**线性趋势**。
- **不适用**：三组无顺序或效应可能呈非单调（如中剂量最优）——请用侧栏「率齐性」。
- **分数**：默认等距 0,1,2；也可改为与剂量成比例的数值，须**严格递增**。
- **样本量**：蒙特卡洛模拟二项数据 + 双侧渐近 p 值；大样本可与统计师用正态近似交叉验证。
"""
                    )
        else:
            if test_type == "精度分析":
                if group_type == "独立两组 (Parallel)":
                    c_p1, c_p2 = st.columns(2)
                    with c_p1:
                        p_treatment_percent = st.number_input("试验组预期有效率 (%)", value=50.0, step=0.1)
                    with c_p2:
                        p_control_percent = st.number_input("对照组预期有效率 (%)", value=40.0, step=0.1)
                    half_width_percent = st.number_input(
                        "率差容许误差半宽 w (%)",
                        value=10.0,
                        step=0.5,
                        help="指两组率差 (pT−pC) 的置信区间半宽，例如 10 表示半宽 0.10",
                    )
                else:
                    p_treatment_percent = st.number_input("预期有效率 P (%)", value=40.0, step=0.1)
                    p_control_percent = 0.0
                    half_width_percent = st.number_input(
                        "容许误差半宽 w (%)",
                        value=10.0,
                        step=0.5,
                        help="例如 10 表示 ±10%，即 95%CI 半宽不超过 0.10",
                    )
                margin_percent = 0.0
                direction = "higher"
            else:
                c_p1, c_p2 = st.columns(2)
                with c_p1:
                    p_treatment_percent = st.number_input("试验组预期有效率 (%)", value=60.0, step=0.1)
                label_ctrl = "对照组预期有效率 (%)" if group_type == "独立两组 (Parallel)" else "历史目标有效率 P0 (%)"
                with c_p2:
                    p_control_percent = st.number_input(label_ctrl, value=40.0, step=0.1)

                # 低值更优：例如不良事件发生率、复发率等（越低越好）
                if group_type == "独立两组 (Parallel)" and test_type == "优效性检验":
                    direction_label = st.selectbox(
                        "指标方向（优效）",
                        ["高值更优（有效率更高）", "低值更优（事件率更低）"],
                        index=0,
                    )
                    direction = "lower" if direction_label.startswith("低值更优") else "higher"
                else:
                    direction = "higher"

                if test_type == "差异性检验":
                    margin_percent = 0.0
                else:
                    margin_label = "非劣界值 Δ" if test_type == "非劣效性检验" else "等效/优效界值 Δ"
                    margin_percent = st.number_input(f"{margin_label} (百分点)", value=5.0, step=0.1)

                if group_type == "独立两组 (Parallel)":
                    st.caption(
                        "双臂率：差异性在期望频数≥5时常用合并方差 Z（与 Pearson χ² 等价）；"
                        "极端率/小样本倾向 Fisher 或反正弦。优/非劣/等效以 FM Score 为药政常用基准，"
                        "等效性对应 FM TOST；可选 Wald 简化或 Chan 类无条件（近似）。精度分析见 Wald / Newcombe。"
                    )
                    with st.expander("双臂率算法说明（ICH / 临床研究常用）", expanded=False):
                        st.markdown(
                            """
| 检验目标 | 场景 | 推荐算法 | 备选/补充 |
|---------|------|----------|-----------|
| 差异性（Δ=0） | 期望频数 ≥5（Cochran） | **Pooled Z（合并方差）** | 与 Pearson χ² 等价 |
|  | 极端率 / 小样本 | **Fisher 精确** | 控制 I 类错误 |
| 优效 / 非劣 / 等效 | 期望频数较充分 | **FM Score（边界投影）** | 药政申报常用 |
|  | 极端率 / 极小样本 | **Chan 无条件精确**（本工具为 FM 近似） | StatXact 等复核 |
|  | 大样本简化 | **Wald（非合并边界方差）** | 与 FM 接近 |
| 等效性（TOST） | 期望频数较充分 | **FM TOST** | 两个单侧 FM |
| 精度分析（CI） | 率差置信区间 | **Wald** / **Newcombe-Wilson** | 极端率更稳 |

* **自动（差异性）**：按合并方差估计样本量后，若四格最小期望频数 ≥5（Cochran）→ Pooled Z；否则 → Fisher 精确。罕见病若 n 足够大仍可走 Pooled Z；极端率且期望频数不足时请手动选 Fisher 或反正弦。*
"""
                        )

        st.markdown('<div class="section-title">统计学设计</div>', unsafe_allow_html=True)
        c_stat1, c_stat2 = st.columns(2)
        with c_stat1:
            if test_type == "精度分析":
                alpha = st.selectbox("显著性水平 α（置信度 = 1−α）", [0.05, 0.025, 0.01, 0.10], index=0)
                power_percent = 80.0
            else:
                power_percent = st.number_input("检验效能 (Power, %)", value=80.0, step=1.0)
                alpha = st.selectbox("显著性水平 (Alpha)", [0.05, 0.025, 0.01], index=0)
            if group_type == "独立两组 (Parallel)":
                alloc_ratio_str = st.selectbox("分配比例 (试验:对照)", ["1:1", "2:1", "3:1", "1:2"], index=0)
            else:
                alloc_ratio_str = "1:1"
        with c_stat2:
            dropout_percent = st.number_input("预计脱落率 (%)", value=20.0, step=1.0)

            if test_type == "精度分析":
                if group_type == "独立两组 (Parallel)":
                    precision_method = st.selectbox(
                        "计算方法（率差 CI）",
                        ["Wald 近似（简单直观）", "Newcombe-Wilson Score（极端率更稳）"],
                        index=0,
                    )
                else:
                    precision_method = st.selectbox(
                        "计算方法",
                        ["正态近似法 (Wald)", "确切概率法 (Clopper-Pearson)"],
                        index=0,
                    )
                alpha_side_label = "双侧"
            elif variable_type == "分类变量 (率)" and group_type == "独立两组 (Parallel)" and test_type != "精度分析":
                if test_type == "差异性检验":
                    _diff_lbl = st.selectbox(
                        "计算方法（差异性）",
                        [
                            "自动（推荐：期望频数≥5→合并方差 Z，否则→Fisher 精确）",
                            "合并方差 Z（与 Pearson 卡方等价）",
                            "Fisher 精确检验",
                            "反正弦变换（极端率备选）",
                        ],
                        index=0,
                        help="自动：按合并方差估计的样本量下，若四格最小期望频数≥5 用 Pooled Z，否则 Fisher。"
                        " 极端率可选手动反正弦或 Fisher。",
                    )
                    _dmap = {
                        "自动（推荐：期望频数≥5→合并方差 Z，否则→Fisher 精确）": None,
                        "合并方差 Z（与 Pearson 卡方等价）": "pooled_z",
                        "Fisher 精确检验": "fisher_exact",
                        "反正弦变换（极端率备选）": "arcsine",
                    }
                    two_arm_diff_method = _dmap[_diff_lbl]
                else:
                    _marg_lbl = st.selectbox(
                        "计算方法（优效/非劣/等效）",
                        [
                            "自动（FM Score）",
                            "FM Score（Farrington–Manning，药政常用）",
                            "Wald 非合并（大样本简化）",
                            "Chan 无条件精确（FM 近似，需复核）",
                        ],
                        index=0,
                    )
                    _mmap = {
                        "自动（FM Score）": None,
                        "FM Score（Farrington–Manning，药政常用）": "fm_score",
                        "Wald 非合并（大样本简化）": "wald_unpooled",
                        "Chan 无条件精确（FM 近似，需复核）": "unconditional_chan",
                    }
                    two_arm_margin_method = _mmap[_marg_lbl]
            elif variable_type == "分类变量 (率)" and group_type == "单臂设计 (Single Arm)":
                method_label = st.selectbox("计算方法", ["正态近似 z 检验（Score）", "精确二项检验（Exact）"], index=0)
                single_arm_method = "exact" if method_label.startswith("精确") else "score"
            elif variable_type == "分类变量 (率)" and group_type == "配对两组 (Paired)":
                if test_type == "差异性检验":
                    _pr_lbl = st.selectbox(
                        "计算方法（配对率）",
                        [
                            "自动（正常→McNemar，小样本/discordant 少→Exact）",
                            "McNemar 检验（正态近似）",
                            "Exact McNemar（小样本 / discordant 很少）",
                        ],
                        index=0,
                    )
                    _pr_map = {
                        "自动（正常→McNemar，小样本/discordant 少→Exact）": None,
                        "McNemar 检验（正态近似）": "nam_score",
                        "Exact McNemar（小样本 / discordant 很少）": "exact_mcnemar",
                    }
                elif test_type == "等效性检验":
                    _pr_lbl = st.selectbox(
                        "计算方法（配对率）",
                        [
                            "Nam score TOST（推荐）",
                            "paired Wald-type（快速近似，非主方法）",
                        ],
                        index=0,
                    )
                    _pr_map = {
                        "Nam score TOST（推荐）": "nam_score",
                        "paired Wald-type（快速近似，非主方法）": "wald_blackwelder",
                    }
                else:
                    _pr_lbl = st.selectbox(
                        "计算方法（配对率）",
                        [
                            "Nam score 方法（推荐）",
                            "paired Wald-type（快速近似，非主方法）",
                            "exact conditional paired（极小样本/敏感性）",
                        ],
                        index=0,
                    )
                    _pr_map = {
                        "Nam score 方法（推荐）": "nam_score",
                        "paired Wald-type（快速近似，非主方法）": "wald_blackwelder",
                        "exact conditional paired（极小样本/敏感性）": "exact_conditional_paired",
                    }
                paired_rate_method = _pr_map[_pr_lbl]
            elif (
                variable_type == "分类变量 (率)"
                and group_type == "三组设计 (Three Arm)"
                and test_type == "差异性检验（率齐性）"
            ):
                _ta_lbl = st.selectbox(
                    "计算方法（三组率齐性）",
                    [
                        "自动（推荐：期望频数充分→Pearson χ²，否则→蒙特卡洛）",
                        "Pearson χ² 齐性（非中心 χ²，常规）",
                        "蒙特卡洛模拟（小样本/极端率/复核）",
                    ],
                    index=0,
                    help="齐性检验；等额三组。蒙特卡洛较慢但适合极小样本或方法学复核。",
                )
                _ta_map = {
                    "自动（推荐：期望频数充分→Pearson χ²，否则→蒙特卡洛）": None,
                    "Pearson χ² 齐性（非中心 χ²，常规）": "pearson",
                    "蒙特卡洛模拟（小样本/极端率/复核）": "monte_carlo",
                }
                three_arm_diff_method = _ta_map[_ta_lbl]

            if test_type != "精度分析":
                if test_type in [
                    "差异性检验",
                    "等效性检验",
                    "差异性检验（率齐性）",
                    "有序趋势检验（Cochran-Armitage）",
                ]:
                    alpha_side_label = "双侧"
                    st.text_input("检验方向", value="双侧检验 (Two-sided)", disabled=True)
                elif variable_type == "生存时间 (HR)" or group_type == "Simon 二阶段" or test_type in ["优效性检验", "非劣效性检验"]:
                    alpha_side_label = "单侧"
                    st.text_input("检验方向", value="单侧检验 (One-sided)", disabled=True)
                else:
                    alpha_side_label = "单侧"
                    st.text_input("检验方向", value="单侧检验 (One-sided)", disabled=True)

    # ------------------ 后端调度区 ------------------
    proportion_calculator = ProportionCalculator()
    continuous_calculator = ContinuousCalculator()
    survival_calculator = SurvivalCalculator()
    error_message, result, simon_res = None, None, None
    alloc_ratio_val = eval(alloc_ratio_str.replace(":", "/"))

    try:
        if variable_type == "分类变量 (率)" and group_type == "Simon 二阶段":
            simon_res = proportion_calculator.compute_simon_two_stage(
                p0=p_control_percent / 100,
                p1=p_treatment_percent / 100,
                alpha=alpha,
                power=power_percent / 100,
                dropout_rate=dropout_percent / 100,
                design_type=simon_design_type or "optimal",
            )
        elif variable_type == "连续变量 (均值)":
            if group_type == "独立两组 (Parallel)":
                result = continuous_calculator.compute_two_arm(
                    mean_t=mean_t,
                    mean_c=mean_c,
                    sd_t=sd_t,
                    sd_c=sd_c,
                    alpha=alpha,
                    power=power_percent / 100,
                    margin=margin_percent,
                    allocation_ratio=alloc_ratio_val,
                    dropout_rate=dropout_percent / 100,
                    test_type=test_type,
                )
            elif group_type == "单臂设计 (Single Arm)":
                result = continuous_calculator.compute_single_arm(
                    mean_t=mean_t,
                    mean_ref=mean_ref,
                    sd=sd_single,
                    alpha=alpha,
                    power=power_percent / 100,
                    margin=margin_percent,
                    dropout_rate=dropout_percent / 100,
                    test_type=test_type,
                    method=_c.get("single_arm_method", "z"),
                )
            elif group_type == "配对两组 (Paired)":
                result = continuous_calculator.compute_paired(
                    mean_diff=mean_diff,
                    sd_diff=sd_diff,
                    alpha=alpha,
                    power=power_percent / 100,
                    margin=margin_percent,
                    dropout_rate=dropout_percent / 100,
                    test_type=test_type,
                )
            else:
                raise ValueError("该连续变量分组设计暂未开放。")
        elif variable_type == "生存时间 (HR)":
            if group_type == "独立两组 (Parallel)":
                result = survival_calculator.compute_two_arm(
                    hr_input_mode="median",
                    hr_direct=None,
                    median_t=median_t,
                    median_c=median_c,
                    alpha=alpha,
                    power=power_percent / 100,
                    allocation_ratio=alloc_ratio_val,
                    dropout_rate=dropout_percent / 100,
                    test_type=test_type,
                    margin_hr=None,
                    equiv_delta=None,
                    event_mode="direct",
                    event_rate=event_rate_percent / 100,
                    T_a=None,
                    T_f=None,
                    lr_method="schoenfeld",
                )
            else:
                result = survival_calculator.compute_single_arm(
                    hr_input_mode="median",
                    hr_direct=None,
                    median_t=median_t,
                    median_ref=median_ref,
                    alpha=alpha,
                    power=power_percent / 100,
                    dropout_rate=dropout_percent / 100,
                    test_type=test_type,
                    margin_hr=None,
                    equiv_delta=None,
                    event_mode="direct",
                    event_rate=event_rate_percent / 100,
                    T_a=None,
                    T_f=None,
                )
        else:
            if (
                variable_type == "分类变量 (率)"
                and group_type == "独立两组 (Parallel)"
                and test_type == "精度分析"
            ):
                pm = "newcombe_wilson" if precision_method.startswith("Newcombe") else "wald"
                prec = proportion_calculator.compute_two_arm_precision(
                    p_t=p_treatment_percent / 100,
                    p_c=p_control_percent / 100,
                    r=alloc_ratio_val,
                    ci_level=1 - alpha,
                    half_width=half_width_percent / 100,
                    method=pm,
                )
                nt = int(prec["n_treatment"])
                nc = int(prec["n_control"])
                ntd = int(math.ceil(nt / (1 - dropout_percent / 100)))
                ncd = int(math.ceil(nc / (1 - dropout_percent / 100)))
                result = {
                    "n_treatment": nt,
                    "n_control": nc,
                    "total_sample_size": nt + nc,
                    "n_treatment_with_dropout": ntd,
                    "n_control_with_dropout": ncd,
                    "total_sample_size_with_dropout": ntd + ncd,
                    "method": prec["method"],
                    "precision": prec,
                }
            elif (
                variable_type == "分类变量 (率)"
                and group_type == "单臂设计 (Single Arm)"
                and test_type == "精度分析"
            ):
                prec_method = "clopper_pearson" if precision_method.startswith("确切") else "wald"
                prec = proportion_calculator.compute_single_arm_precision(
                    p=p_treatment_percent / 100,
                    ci_level=1 - alpha,
                    half_width=half_width_percent / 100,
                    method=prec_method,
                )
                n = int(prec["n"])
                n_drop = int(math.ceil(n / (1 - dropout_percent / 100)))
                result = {
                    "n_treatment": n,
                    "n_control": 0,
                    "total_sample_size": n,
                    "n_treatment_with_dropout": n_drop,
                    "n_control_with_dropout": 0,
                    "total_sample_size_with_dropout": n_drop,
                    "method": prec["method"],
                    "precision": prec,
                }
            elif (
                variable_type == "分类变量 (率)"
                and group_type == "配对两组 (Paired)"
            ):
                result = proportion_calculator.compute_paired_proportion(
                    p_discordant_treat=p10_pct / 100.0,
                    p_discordant_ctrl=p01_pct / 100.0,
                    alpha=alpha,
                    power=power_percent / 100.0,
                    dropout_rate=dropout_percent / 100.0,
                    test_type=test_type,
                    margin=margin_percent / 100.0 if test_type != "差异性检验" else 0.0,
                    method=paired_rate_method,
                )
            elif variable_type == "分类变量 (率)" and group_type == "三组设计 (Three Arm)":
                if test_type == "差异性检验（率齐性）":
                    result = proportion_calculator.compute_three_arm_proportion(
                        p1_pct / 100.0,
                        p2_pct / 100.0,
                        p3_pct / 100.0,
                        alpha=alpha,
                        power=power_percent / 100.0,
                        dropout_rate=dropout_percent / 100.0,
                        method=three_arm_diff_method,
                    )
                else:
                    result = proportion_calculator.compute_three_arm_cochran_armitage(
                        p1_pct / 100.0,
                        p2_pct / 100.0,
                        p3_pct / 100.0,
                        alpha=alpha,
                        power=power_percent / 100.0,
                        dropout_rate=dropout_percent / 100.0,
                        scores=(ca_t1, ca_t2, ca_t3),
                    )
            else:
                params = ProportionInput(
                    p_treatment=p_treatment_percent / 100,
                    p_control=p_control_percent / 100,
                    alpha=alpha,
                    power=power_percent / 100,
                    test_type=test_type,
                    allocation_ratio=alloc_ratio_val,
                    dropout_rate=dropout_percent / 100,
                    margin=margin_percent / 100,
                    direction=direction if group_type == "独立两组 (Parallel)" else "higher",
                    two_arm_diff_method=two_arm_diff_method,
                    two_arm_margin_method=two_arm_margin_method,
                )
                if group_type == "独立两组 (Parallel)":
                    result = proportion_calculator.compute_two_arm(params)
                else:
                    result = proportion_calculator.compute_single_arm(params, method=single_arm_method)
    except Exception as e:
        error_message = str(e)

    # ------------------ 右侧结果报告区 ------------------
    with col_result:
        if error_message:
            st.error(f"计算冲突: {error_message}")
            return

        # 渲染 Simon 二阶段结果
        if group_type == "Simon 二阶段" and simon_res:
            st.markdown(
                f"""
                <div class="result-kpi-card">
                    <div style="font-size: 16px; font-weight: 500; text-align: center; opacity: 0.95;">Simon {simon_res['design_type'].capitalize()} 目标招募样本量</div>
                    <div class="group-split">
                        <div class="group-box">
                            <div class="group-num">{simon_res['n1']}</div>
                            <div class="group-name">第一阶段</div>
                        </div>
                        <div class="group-box">
                            <div class="group-num">{simon_res['n2']}</div>
                            <div class="group-name">第二阶段</div>
                        </div>
                    </div>
                    <div class="total-box">
                        理论总样本量：{simon_res['n']} 例；实际计划入组：{simon_res['total_sample_size_with_dropout']} 例
                        <span style="font-weight: 400; font-size: 15px; opacity: 0.9;">(已计入 {dropout_percent}% 脱落率)</span>
                    </div>
                </div>
                """, unsafe_allow_html=True
            )

            st.caption(
                f"设计参数：r1={simon_res['r1']}，r={simon_res['r']}，"
                f"PET(p0)={simon_res['pet0']*100:.1f}%，EN(p0)={simon_res['en0']:.1f}，"
                f"实际α={simon_res['actual_alpha']:.4f}，实际Power={simon_res['actual_power']:.4f}"
            )
            st.caption(
                f"对照设计：Optimal(n={simon_res['optimal_design']['n']}, EN0={simon_res['optimal_design']['en0']:.1f})；"
                f"Minimax(n={simon_res['minimax_design']['n']}, EN0={simon_res['minimax_design']['en0']:.1f})"
            )
            
            st.markdown('<div class="section-title" style="margin-top:0;">样本量计算方法</div>', unsafe_allow_html=True)
            cn_text = (
                f"本研究采用 Simon 二阶段 {simon_res['design_type']} 设计。设定药物无效界值 P0 为 {p_control_percent}%，"
                f"预期目标有效率 P1 为 {p_treatment_percent}%。在单侧显著性水平 α = {alpha} 且检验效能 1-β = {power_percent}% 的条件下，"
                f"第一阶段计划入组 {simon_res['n1']} 例受试者。若在第一阶段观察到 ≤ {simon_res['r1']} 例有效，则因无效提前终止试验"
                f"（早期停止概率 PET = {simon_res['pet0']*100:.1f}%）；若观察到 > {simon_res['r1']} 例有效，则进入第二阶段，"
                f"继续入组 {simon_res['n2']} 例受试者。<br><br>"
                f"最终，若在总计理论 {simon_res['n']} 例受试者中观察到 > {simon_res['r']} 例有效，则拒绝原假设，"
                f"认为该药物具有临床疗效。考虑 {dropout_percent}% 的脱落率，实际计划入组最大样本量为 {simon_res['total_sample_size_with_dropout']} 例。"
            )
            st.markdown(f"<div class='statement-box'><p>{cn_text}</p></div>", unsafe_allow_html=True)

        # 渲染 单臂/双臂结果
        elif result:
            total_n = result["total_sample_size_with_dropout"]
            n_treat = result.get("n_treatment_with_dropout", 0)
            n_ctrl = result.get("n_control_with_dropout", 0)
            theo_total = result["total_sample_size"]
            theo_treat = int(result.get("n_treatment", theo_total / (1 + 1/alloc_ratio_val)))
            theo_ctrl = int(result.get("n_control", theo_total - theo_treat))

            # 三组率（等额）
            if variable_type == "分类变量 (率)" and group_type == "三组设计 (Three Arm)":
                theo_pg = int(result["n_per_group"])
                npg_drop = int(result["n_per_group_with_dropout"])
                st.markdown(
                    f"""
                    <div class="result-kpi-card">
                        <div style="font-size: 16px; font-weight: 500; text-align: center; opacity: 0.95;">目标招募总样本量</div>
                        <div style="text-align: center; margin: 30px 0;">
                            <span style="font-size: 64px; font-weight: 700;">{total_n}</span>
                        </div>
                        <div class="total-box">理论所需：{theo_total} 例（三组等额，每组 {theo_pg} 例）
                            <span style="font-weight: 400; font-size: 15px; opacity: 0.9;">(已计入 {dropout_percent}% 脱落率)</span>
                        </div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
                st.caption(f"计划入组：每组约 {npg_drop} 例（含脱落）。")
                if result.get("method"):
                    st.caption(result["method"])
                if result.get("rate_note"):
                    st.caption(result["rate_note"])
                if result.get("min_expected_cell_approx") is not None:
                    st.caption(
                        f"按假设率估计的齐性检验最小期望频数（每组阳性/阴性中较小者）≈ {result['min_expected_cell_approx']:.2f}"
                    )
                if result.get("actual_power_mc") is not None:
                    st.caption(
                        f"蒙特卡洛估计功效 ≈ {result['actual_power_mc'] * 100:.2f}%（{result.get('mc_replications', '')} 次重复）"
                    )
                st.markdown('<div class="section-title" style="margin-top:0;">样本量计算方法</div>', unsafe_allow_html=True)
                if test_type == "差异性检验（率齐性）":
                    cn_three_1 = (
                        f"本研究为三组独立样本二分类终点**率齐性**分析。假定三组有效率分别约为 π₁={p1_pct:.1f}%、"
                        f"π₂={p2_pct:.1f}%、π₃={p3_pct:.1f}%。"
                        f"采用双侧检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%，三组等额分配。"
                    )
                else:
                    sc = result.get("ca_scores") or (0.0, 1.0, 2.0)
                    cn_three_1 = (
                        f"本研究为三组有序二分类终点的**线性趋势**分析（Cochran–Armitage，双侧）。"
                        f"假定三组有效率分别约为 π₁={p1_pct:.1f}%、π₂={p2_pct:.1f}%、π₃={p3_pct:.1f}%，"
                        f"对应趋势分数为 ({sc[0]:g}, {sc[1]:g}, {sc[2]:g})。"
                        f"采用双侧检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%，三组等额分配。"
                    )
                cn_three_2 = (
                    f"样本量采用【{result.get('method', '')}】估算。"
                    f"理论所需总样本量为 {theo_total} 例（每组 {theo_pg} 例）；考虑约 {dropout_percent:.1f}% 脱落率，计划总入组 {total_n} 例。"
                )
                st.markdown(
                    f"<div class='statement-box'><p>{cn_three_1}</p><p>{cn_three_2}</p></div>",
                    unsafe_allow_html=True,
                )
            # 单臂UI
            elif group_type in ["单臂设计 (Single Arm)", "配对两组 (Paired)"]:
                _kpi_title = (
                    "目标配对受试者数"
                    if variable_type == "分类变量 (率)" and group_type == "配对两组 (Paired)"
                    else "目标招募样本量"
                )
                _kpi_sub = (
                    "每人提供配对二分类（两次测量）"
                    if variable_type == "分类变量 (率)" and group_type == "配对两组 (Paired)"
                    else ""
                )
                st.markdown(
                    f"""
                    <div class="result-kpi-card">
                        <div style="font-size: 16px; font-weight: 500; text-align: center; opacity: 0.95;">{_kpi_title}</div>
                        <div style="text-align: center; margin: 30px 0;">
                            <span style="font-size: 64px; font-weight: 700;">{total_n}</span>
                        </div>
                        <div class="total-box">理论所需：{theo_total} 例 <span style="font-weight: 400; font-size: 15px; opacity: 0.9;">(已计入 {dropout_percent}% 脱落率)</span></div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
                if _kpi_sub:
                    st.caption(_kpi_sub)

                st.markdown('<div class="section-title" style="margin-top:0;">样本量计算方法</div>', unsafe_allow_html=True)
                if variable_type == "分类变量 (率)":
                    if (
                        group_type == "配对两组 (Paired)"
                        and test_type != "精度分析"
                    ):
                        margin_txt = (
                            f"，设定率差界值 Δ = {margin_percent:.1f}%"
                            if test_type != "差异性检验"
                            else ""
                        )
                        cn_text_p1 = (
                            f"本研究为配对两组二分类终点{test_type[:2]}设计。假定不一致配对比例 π₁₀≈{p10_pct:.2f}%、"
                            f"π₀₁≈{p01_pct:.2f}%（ψ=π₁₀−π₀₁）{margin_txt}。"
                            f"采用{alpha_side_label}检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                        )
                        cn_text_p2 = (
                            f"样本量采用【{result.get('method', '')}】计算。"
                            f"理论所需配对受试者数为 {theo_total} 例；考虑约 {dropout_percent:.1f}% 脱落率，计划入组 {total_n} 例。"
                        )
                        if isinstance(result, dict) and result.get("rate_note"):
                            st.caption(result["rate_note"])
                    elif test_type == "精度分析":
                        cn_text_p1 = (
                            f"本研究为单臂率精度分析。假定目标率约为 {p_treatment_percent:.1f}%，"
                            f"设定 {(1 - alpha) * 100:.0f}% 置信区间半宽不超过 {half_width_percent:.1f}%。"
                        )
                        cn_text_p2 = (
                            f"样本量采用【{result.get('method', '')}】计算。"
                            f"理论最小样本量为 {theo_total} 例；考虑约 {dropout_percent:.1f}% 脱落率，计划入组 {total_n} 例。"
                        )
                    else:
                        margin_txt = f"，设定{test_type[:2]}界值 Δ = {margin_percent:.1f}%" if test_type != "差异性检验" else ""
                        method_line = result.get("method", "单臂率：Score 正态近似") if isinstance(result, dict) else "单臂率：Score 正态近似"
                        cn_text_p1 = (
                            f"本研究为单臂{test_type[:2]}设计。根据临床预期，假定试验组有效率为 {p_treatment_percent:.1f}%，"
                            f"历史对照目标值为 {p_control_percent:.1f}%{margin_txt}。采用{alpha_side_label}检验，"
                            f"显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                        )
                        cn_text_p2 = (
                            f"样本量采用【{method_line}】进行估计。在上述参数条件下，理论所需总样本量为 {theo_total} 例。"
                            f"考虑约 {dropout_percent:.1f}% 的脱落率，为保证最终分析所需样本量，实际计划总入组 {total_n} 例。"
                        )
                    # 若为精确二项，额外展示实际 α / power 与搜索起点信息，证明确实跑了 Exact
                    if isinstance(result, dict) and "actual_alpha_at_p0" in result:
                        st.caption(
                            f"Exact 校核：实际α(p0)≈{result['actual_alpha_at_p0']:.4f}，"
                            f"实际Power(p1)≈{result.get('actual_power_at_p1', float('nan')):.4f}；"
                            f"搜索起点 n0={result.get('n0_start','-')}，n_start={result.get('n_start','-')}。"
                        )
                elif variable_type == "连续变量 (均值)":
                    if group_type == "配对两组 (Paired)":
                        margin_txt = f"，界值 Δ = {margin_percent:.2f}" if test_type != "差异性检验" else ""
                        cn_text_p1 = (
                            f"本研究为配对两组连续变量{test_type[:2]}设计。假定配对差值均值为 {mean_diff:.2f}，"
                            f"配对差值标准差 SD_diff 为 {sd_diff:.2f}{margin_txt}。采用{alpha_side_label}检验，"
                            f"显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                        )
                        cn_text_p2 = (
                            f"样本量基于配对差值的正态近似方法进行计算。在上述参数条件下，理论所需样本量为 {theo_total} 例。"
                            f"考虑约 {dropout_percent:.1f}% 的脱落率，为保证最终分析所需样本量，实际计划入组 {total_n} 例。"
                        )
                    else:
                        margin_txt = f"，界值 Δ = {margin_percent:.2f}" if test_type != "差异性检验" else ""
                        cn_text_p1 = (
                            f"本研究为单臂连续变量{test_type[:2]}设计。假定试验组均值为 {mean_t:.2f}，"
                            f"参考均值为 {mean_ref:.2f}，标准差 SD 为 {sd_single:.2f}{margin_txt}。"
                            f"采用{alpha_side_label}检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                        )
                        _meth = (
                            result.get("method", "单样本均值：正态近似 z")
                            if isinstance(result, dict)
                            else "单样本均值：正态近似 z"
                        )
                        cn_text_p2 = (
                            f"样本量采用【{_meth}】。在上述参数条件下，理论所需样本量为 {theo_total} 例。"
                            f"考虑约 {dropout_percent:.1f}% 的脱落率，为保证最终分析所需样本量，实际计划入组 {total_n} 例。"
                        )
                else:
                    hr_show = result.get("hr", None)
                    event_show = result.get("events", None)
                    hr_clause = f"，对应预期 HR≈{hr_show:.3f}" if hr_show is not None else ""
                    evt_clause = f"预计所需事件数约 {event_show} 例，" if event_show is not None else ""
                    cn_text_p1 = (
                        f"本研究为单臂生存终点{test_type[:2]}设计。假定试验组中位生存期 {median_t:.1f} 月，"
                        f"参考中位生存期 {median_ref:.1f} 月{hr_clause}。采用{alpha_side_label}检验，"
                        f"显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                    )
                    cn_text_p2 = (
                        f"{evt_clause}预计事件发生率为 {event_rate_percent:.1f}% 时，理论所需样本量为 {theo_total} 例。"
                        f"考虑约 {dropout_percent:.1f}% 的脱落率，实际计划入组 {total_n} 例。"
                    )
                st.markdown(f"<div class='statement-box'><p>{cn_text_p1}</p><p>{cn_text_p2}</p></div>", unsafe_allow_html=True)

            # 双臂UI
            else:
                st.markdown(
                    f"""
                    <div class="result-kpi-card">
                        <div style="font-size: 16px; font-weight: 500; text-align: center; opacity: 0.95; letter-spacing: 0.5px;">目标招募样本量</div>
                        <div class="group-split">
                            <div class="group-box">
                                <div class="group-num">{n_treat}</div>
                                <div class="group-name">试验组</div>
                            </div>
                            <div class="group-box">
                                <div class="group-num">{n_ctrl}</div>
                                <div class="group-name">对照组</div>
                            </div>
                        </div>
                        <div class="total-box">
                            总计招募：{total_n} 例 <span style="font-weight: 400; font-size: 15px; opacity: 0.9;">(已计入 {dropout_percent}% 脱落率)</span>
                        </div>
                    </div>
                    """, unsafe_allow_html=True
                )

                if variable_type == "分类变量 (率)" and group_type == "独立两组 (Parallel)" and result:
                    if result.get("method"):
                        st.caption(result["method"])
                    if result.get("expected_min_cell") is not None:
                        st.caption(
                            f"基于假设率估计的四格最小期望频数 ≈ {result['expected_min_cell']:.2f}"
                        )
                    if result.get("rate_note"):
                        st.caption(result["rate_note"])
                
                st.markdown('<div class="section-title" style="margin-top:0;">样本量计算方法</div>', unsafe_allow_html=True)
                alloc_txt = f"即每组 {theo_treat} 例" if alloc_ratio_val == 1.0 else f"即试验组 {theo_treat} 例，对照组 {theo_ctrl} 例"

                if variable_type == "分类变量 (率)":
                    if test_type == "精度分析":
                        cn_text_p1 = (
                            f"本研究为两臂率差精度分析。假定试验组与对照组有效率分别约为 {p_treatment_percent:.1f}%、{p_control_percent:.1f}%，"
                            f"设定 {(1 - alpha) * 100:.0f}% 置信水平下率差置信区间半宽不超过 {half_width_percent:.1f}%。"
                        )
                        cn_text_p2 = (
                            f"样本量采用【{result.get('method', '')}】计算。"
                            f"理论最小总样本量为 {theo_total} 例（{alloc_txt}）；考虑约 {dropout_percent:.1f}% 脱落率，计划入组 {total_n} 例。"
                        )
                    else:
                        margin_txt = f"，设定{test_type[:2]}界值 Δ = {margin_percent:.1f}%" if test_type != "差异性检验" else ""
                        if result and result.get("rate_zone") == "fisher":
                            meth_line = "Fisher 精确检验（两独立样本率，双侧）"
                        elif result and result.get("rate_zone") == "extreme" and test_type == "差异性检验":
                            meth_line = "反正弦变换（方差稳定化）法（两独立样本率差）"
                        elif (
                            result
                            and result.get("rate_zone") == "extreme_pooled"
                            and test_type == "差异性检验"
                        ):
                            meth_line = "合并方差 Z（与 Pearson χ² 等价；率偏斜但期望频数仍 ≥5）"
                        elif (
                            result
                            and result.get("rate_zone") == "extreme_pooled_sparse"
                            and test_type == "差异性检验"
                        ):
                            meth_line = "合并方差 Z；期望频数不足 5，建议分析阶段 Fisher 或改用自动/反正弦"
                        elif result and result.get("rate_zone") == "extreme" and test_type in [
                            "优效性检验",
                            "非劣效性检验",
                            "等效性检验",
                        ]:
                            meth_line = "Farrington–Manning score / FM Score（两独立样本率差，H₀ 边界投影）"
                        elif result and result.get("rate_zone") == "wald":
                            meth_line = "Wald 型非合并方差（H₀ 边界直接代入）"
                        elif test_type == "等效性检验":
                            meth_line = "FM TOST（两个单侧 FM Score，正态近似）"
                        else:
                            meth_line = "正态近似法（两独立样本率差）"

                        cn_text_p1 = (
                            f"主要终点为二分类结局（有效率）。本研究采用独立两组{test_type[:2]}设计。"
                            f"根据既往研究及临床预期，试验组有效率预计为 {p_treatment_percent:.1f}%，对照组为 {p_control_percent:.1f}%{margin_txt}。"
                            f"采用{alpha_side_label}检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%，按 {alloc_ratio_str} 比例分配。"
                        )
                        cn_text_p2 = (
                            f"样本量采用【{meth_line}】计算。按上述假设与参数，理论所需总样本量为 {theo_total} 例（{alloc_txt}）。"
                            f"考虑约 {dropout_percent:.1f}% 脱落率，为保证达到理论分析样本量，计划最多入组 {total_n} 例，"
                            f"其中试验组 {n_treat} 例、对照组 {n_ctrl} 例。"
                        )
                elif variable_type == "连续变量 (均值)":
                    margin_txt = f"，界值 Δ = {margin_percent:.2f}" if test_type != "差异性检验" else ""
                    cn_text_p1 = (
                        f"本研究为连续变量终点的独立两组{test_type[:2]}设计。假定试验组均值为 {mean_t:.2f}，"
                        f"对照组均值为 {mean_c:.2f}，试验组标准差 {sd_t:.2f}、对照组标准差 {sd_c:.2f}{margin_txt}。"
                        f"采用{alpha_side_label}检验，显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                    )
                    cn_text_p2 = (
                        f"样本量基于两独立样本均值差的正态近似方法进行计算。在上述参数条件下，理论总样本量为 {theo_total} 例，"
                        f"{alloc_txt}。考虑约 {dropout_percent:.1f}% 的脱落率，实际计划总入组 {total_n} 例，"
                        f"其中试验组 {n_treat} 例，对照组 {n_ctrl} 例。"
                    )
                else:
                    hr_show = result.get("hr", None)
                    event_show = result.get("events", None)
                    hr_clause = f"，对应预期 HR≈{hr_show:.3f}" if hr_show is not None else ""
                    evt_clause = f"所需事件数约 {event_show} 例，" if event_show is not None else ""
                    cn_text_p1 = (
                        f"本研究为生存终点的独立两组{test_type[:2]}设计。假定试验组中位生存期 {median_t:.1f} 月，"
                        f"对照组中位生存期 {median_c:.1f} 月{hr_clause}。采用{alpha_side_label}检验，"
                        f"显著性水平 α = {alpha}，检验效能（1−β）= {power_percent:.1f}%。"
                    )
                    cn_text_p2 = (
                        f"样本量基于 log-rank 近似方法计算，预计事件发生率 {event_rate_percent:.1f}% 时，"
                        f"{evt_clause}理论总样本量为 {theo_total} 例，{alloc_txt}。考虑约 {dropout_percent:.1f}% 脱落率，"
                        f"实际计划总入组 {total_n} 例，其中试验组 {n_treat} 例，对照组 {n_ctrl} 例。"
                    )
                st.markdown(f"<div class='statement-box'><p>{cn_text_p1}</p><p>{cn_text_p2}</p></div>", unsafe_allow_html=True)

if __name__ == "__main__":
    v, g, t, sd = sidebar_layout()
    main_layout(v, g, t, sd)