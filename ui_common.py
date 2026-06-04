"""共用 UI：全局 CSS、可复用的 HTML 片段。"""

import streamlit as st


def inject_morandi_ui() -> None:
    st.markdown(
        """
        <style>
        /* 全局背景 */
        body, .stApp {
            background-color: #F4F2ED;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            color: #3C3C3C;
        }

        #MainMenu, footer {visibility: hidden;}

        .main-title {
            font-size: 32px;
            font-weight: 700;
            color: #3C3C3C;
            margin-bottom: 40px;
            margin-top: -20px;
            letter-spacing: 1px;
            text-align: center;
        }

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

        .result-kpi-card {
            background-color: #788A9A;
            border-radius: 12px;
            padding: 28px 24px;
            color: #FFFFFF;
            box-shadow: 0 8px 16px rgba(120, 138, 154, 0.2);
            margin-bottom: 24px;
        }

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

        .total-box {
            text-align: center;
            font-size: 18px;
            font-weight: 600;
            border-top: 1px dashed rgba(255, 255, 255, 0.4);
            padding-top: 20px;
            letter-spacing: 0.5px;
        }

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

        div[data-testid="stNumberInput"] label p, div[data-testid="stSelectbox"] label p {
            color: #5A5A5A;
            font-size: 13px;
            font-weight: 500;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )



def render_statistics_design(variable_type: str, group_type: str, test_type: str) -> dict:
    """统计学设计区：α、Power、脱落率、分配、各算法下拉、检验方向。"""
    single_arm_method = "score"
    half_width_percent = 10.0
    precision_method = "正态近似法 (Wald)"
    two_arm_diff_method = None
    two_arm_margin_method = None
    paired_rate_method = None
    three_arm_diff_method = None

    st.markdown('<div class="section-title">统计学设计</div>', unsafe_allow_html=True)
    c_stat1, c_stat2 = st.columns(2)
    with c_stat1:
        if test_type == "精度分析":
            alpha = st.slider("显著性水平 α（置信度 = 1−α）", min_value=0.0, max_value=0.2, value=0.05, step=0.001, format="%.3f")
            power_percent = 80.0
        else:
            power_percent = st.number_input("检验效能 (Power, %)", value=80.0, step=1.0)
            alpha = st.slider("显著性水平 (Alpha)", min_value=0.0, max_value=0.2, value=0.05, step=0.001, format="%.3f")
        if group_type == "独立两组 (Parallel)":
            alloc_ratio_str = st.selectbox("分配比例 (试验:对照)", ["1:1", "2:1", "3:1", "1:2"], index=0)
        else:
            alloc_ratio_str = "1:1"
    with c_stat2:
        dropout_percent = st.number_input("预计脱落率 (%)", value=20.0, step=1.0)

        if test_type == "精度分析":
            if variable_type == "连续变量 (均值)" and group_type == "单臂设计 (Single Arm)":
                st.caption("均值精度样本量：双侧 **t 分布半宽法**（t_{α/2, n−1}·SD/√n ≤ 半宽）。")
                precision_method = "continuous_t_halfwidth"
            elif variable_type == "连续变量 (均值)" and group_type == "独立两组 (Parallel)":
                st.caption(
                    "两独立样本均值**差**精度：双侧 **t 分布半宽法**（Welch SE；t_{α/2, df}·SE ≤ 半宽）。"
                )
                precision_method = "continuous_two_arm_t_halfwidth"
            elif variable_type == "连续变量 (均值)" and group_type == "配对两组 (Paired)":
                st.caption(
                    "配对差值**均值**精度：双侧 **t 分布半宽法**（t_{α/2, n−1}·SD_diff/√n；n 为对子数）。"
                )
                precision_method = "continuous_paired_t_halfwidth"
            elif group_type == "独立两组 (Parallel)":
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
                        "合并方差正态近似（与 Pearson χ² 等价，药政常用）",
                        "非合并方差简化公式（国内方案常用）",
                        "Yates 连续性校正 χ²（Fleiss 公式，保守）",
                        "Fisher 确切概率法",
                    ],
                    index=0,
                )
                _dmap = {
                    "合并方差正态近似（与 Pearson χ² 等价，药政常用）": "pooled_z",
                    "非合并方差简化公式（国内方案常用）": "unpooled_z",
                    "Yates 连续性校正 χ²（Fleiss 公式，保守）": "yates_correction",
                    "Fisher 确切概率法": "fisher_exact",
                }
                two_arm_diff_method = _dmap[_diff_lbl]
            else:
                _marg_lbl = st.selectbox(
                    "计算方法（优效/非劣/等效）",
                    [
                        "Farrington–Manning Score（restricted MLE，金标准）",
                        "Wald 正态近似（大样本简化）",
                    ],
                    index=0,
                    help="FM Score 使用约束极大似然估计求 H₀ 边界方差，与 SAS/PASS/East 一致；Wald 为简化法，略保守。",
                )
                _mmap = {
                    "Farrington–Manning Score（restricted MLE，金标准）": "fm_score",
                    "Wald 正态近似（大样本简化）": "wald",
                }
                two_arm_margin_method = _mmap[_marg_lbl]
        elif variable_type == "分类变量 (率)" and group_type == "单臂设计 (Single Arm)":
            method_label = st.selectbox("计算方法", ["正态近似 z 检验（Score）", "精确二项检验（Exact）"], index=0)
            single_arm_method = "exact" if method_label.startswith("精确") else "score"
        elif variable_type == "分类变量 (率)" and group_type == "配对两组 (Paired)":
            if test_type == "差异性检验":
                _pr_lbl = st.selectbox(
                    "计算方法（配对率差异性）",
                    [
                        "McNemar 检验（Nam 正态近似）",
                        "Exact McNemar（小样本 / discordant 很少）",
                    ],
                    index=0,
                )
                _pr_map = {
                    "McNemar 检验（Nam 正态近似）": "nam_score",
                    "Exact McNemar（小样本 / discordant 很少）": "exact_mcnemar",
                }
            elif test_type == "等效性检验":
                _pr_lbl = st.selectbox(
                    "计算方法（配对率等效性）",
                    [
                        "Nam Score TOST（推荐）",
                    ],
                    index=0,
                )
                _pr_map = {
                    "Nam Score TOST（推荐）": "nam_score",
                }
            else:
                _pr_lbl = st.selectbox(
                    "计算方法（配对率优效/非劣效）",
                    [
                        "Nam Score 方法（推荐）",
                        "Exact Conditional Paired（极小样本/敏感性）",
                    ],
                    index=0,
                )
                _pr_map = {
                    "Nam Score 方法（推荐）": "nam_score",
                    "Exact Conditional Paired（极小样本/敏感性）": "exact_conditional_paired",
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
                    "Pearson χ² 齐性（非中心 χ² 法）",
                ],
                index=0,
                help="三组等额分配，基于渐近 Pearson χ² 的非中心参数法。",
            )
            _ta_map = {
                "Pearson χ² 齐性（非中心 χ² 法）": "pearson",
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
            elif variable_type == "生存时间 (HR)" or group_type == "Simon 二阶段" or test_type in [
                "优效性检验",
                "非劣效性检验",
            ]:
                alpha_side_label = "单侧"
                st.text_input("检验方向", value="单侧检验 (One-sided)", disabled=True)
            else:
                alpha_side_label = "单侧"
                st.text_input("检验方向", value="单侧检验 (One-sided)", disabled=True)
        else:
            alpha_side_label = "双侧"

    return {
        "alpha": alpha,
        "power_percent": power_percent,
        "dropout_percent": dropout_percent,
        "alloc_ratio_str": alloc_ratio_str,
        "precision_method": precision_method,
        "single_arm_method": single_arm_method,
        "half_width_percent": half_width_percent,
        "two_arm_diff_method": two_arm_diff_method,
        "two_arm_margin_method": two_arm_margin_method,
        "paired_rate_method": paired_rate_method,
        "three_arm_diff_method": three_arm_diff_method,
        "alpha_side_label": alpha_side_label,
    }


def sidebar_layout():
    """侧边栏：主要终点变量、分组设计、假设检验目标；返回 (variable_type, group_type, test_type, simon_design_type)。"""
    with st.sidebar:
        st.markdown(
            "<h3 style='color: #4A4A4A; font-size: 17px; margin-bottom: 20px;'>试验设计导航</h3>",
            unsafe_allow_html=True,
        )
        variable_type = st.selectbox(
            "主要终点变量",
            ["分类变量 (率)", "连续变量 (均值)", "生存时间 (HR)"],
            index=0,
        )
        simon_design_type = None

        if variable_type == "分类变量 (率)":
            group_type = st.selectbox(
                "试验分组设计",
                [
                    "独立两组 (Parallel)",
                    "单臂设计 (Single Arm)",
                    "配对两组 (Paired)",
                    "三组设计 (Three Arm)",
                    "Simon 二阶段",
                ],
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
                        [
                            "优效性检验",
                            "非劣效性检验",
                            "差异性检验",
                            "等效性检验",
                            "精度分析",
                        ],
                        index=0,
                    )
                else:
                    test_type = st.selectbox(
                        "假设检验目标",
                        [
                            "优效性检验",
                            "非劣效性检验",
                            "差异性检验",
                            "等效性检验",
                        ],
                        index=0,
                    )
        elif variable_type == "连续变量 (均值)":
            group_type = st.selectbox(
                "试验分组设计",
                [
                    "独立两组 (Parallel)",
                    "单臂设计 (Single Arm)",
                    "配对两组 (Paired)",
                    "三组设计 (Three Arm)",
                ],
                index=0,
            )
            if group_type == "三组设计 (Three Arm)":
                test_type = st.selectbox(
                    "假设检验目标",
                    ["差异性检验"],
                    index=0,
                    help="三组连续变量：差异性方法在临床预期中选择（ANOVA、Kruskal-Wallis、有序线性趋势、Dunnett 等）。",
                )
            else:
                if group_type == "单臂设计 (Single Arm)":
                    test_type = st.selectbox(
                        "假设检验目标",
                        [
                            "优效性检验",
                            "非劣效性检验",
                            "差异性检验",
                            "等效性检验",
                            "精度分析",
                        ],
                        index=0,
                    )
                elif group_type == "独立两组 (Parallel)":
                    test_type = st.selectbox(
                        "假设检验目标",
                        [
                            "优效性检验",
                            "非劣效性检验",
                            "差异性检验",
                            "等效性检验",
                            "精度分析",
                        ],
                        index=0,
                    )
                elif group_type == "配对两组 (Paired)":
                    test_type = st.selectbox(
                        "假设检验目标",
                        [
                            "优效性检验",
                            "非劣效性检验",
                            "差异性检验",
                            "等效性检验",
                            "精度分析",
                        ],
                        index=0,
                    )
                else:
                    test_type = st.selectbox(
                        "假设检验目标",
                        ["优效性检验", "非劣效性检验", "差异性检验", "等效性检验"],
                        index=0,
                    )
        else:
            group_type = st.selectbox(
                "试验分组设计",
                ["独立两组 (Parallel)", "单臂设计 (Single Arm)"],
                index=0,
            )
            test_type = st.selectbox(
                "假设检验目标",
                [
                    "优效性检验",
                    "非劣效性检验",
                    "差异性检验",
                    "等效性检验",
                ],
                index=0,
            )

    return variable_type, group_type, test_type, simon_design_type
