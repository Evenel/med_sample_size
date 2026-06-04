# calc_continuous.py
from __future__ import annotations

import math
from scipy.stats import norm, f as f_dist, ncf, nct as nct_dist, t as t_dist

def _welch_df(v_t: float, v_c: float, nt: float, nc: float) -> float:
    s_e2 = v_t / nt + v_c / nc
    if s_e2 <= 1e-15:
        return 1.0
    return (s_e2 ** 2) / ((v_t / nt) ** 2 / max(nt - 1, 1e-9) + (v_c / nc) ** 2 / max(nc - 1, 1e-9))


class ContinuousCalculator:
    def compute_two_arm_precision(
        self,
        sd_t: float,
        sd_c: float,
        half_width: float,
        alpha: float,
        allocation_ratio: float,
        dropout_rate: float,
    ) -> dict:
        """两独立样本均值差精度：双侧 CI 半宽 ≤ w；SE 用 Welch，临界值用 Welch-Satterthwaite df 下的 t。"""
        if sd_t <= 0 or sd_c <= 0 or half_width <= 0:
            raise ValueError("标准差与均值差 CI 半宽须为正。")
        r = allocation_ratio
        if r <= 0:
            raise ValueError("分配比例须为正。")
        nt = max(2.0, (norm.ppf(1 - alpha / 2) * math.sqrt(sd_t**2 + sd_c**2 / r) / half_width) ** 2)
        nt = max(2, int(math.ceil(nt)))
        for _ in range(500000):
            nc = max(2, int(math.ceil(nt * r)))
            se = math.sqrt(sd_t**2 / nt + sd_c**2 / nc)
            df = _welch_df(sd_t**2, sd_c**2, float(nt), float(nc))
            t_crit = t_dist.ppf(1 - alpha / 2, df)
            if t_crit * se <= half_width:
                break
            nt += 1
        else:
            raise ValueError("在合理范围内未能找到满足精度要求的样本量。")

        nt_ceil, nc_ceil = nt, nc
        nt_drop = int(math.ceil(nt_ceil / (1 - dropout_rate)))
        nc_drop = int(math.ceil(nc_ceil / (1 - dropout_rate)))
        return {
            "n_treatment": nt_ceil,
            "n_control": nc_ceil,
            "total_sample_size": nt_ceil + nc_ceil,
            "n_treatment_with_dropout": nt_drop,
            "n_control_with_dropout": nc_drop,
            "total_sample_size_with_dropout": nt_drop + nc_drop,
            "method": "两独立样本均值差精度：t 分布半宽法（Welch SE，双侧 CI）",
        }

    def compute_two_arm(
        self,
        mean_t: float,
        mean_c: float,
        sd_t: float,
        sd_c: float,
        alpha: float,
        power: float,
        margin: float,
        allocation_ratio: float,
        dropout_rate: float,
        test_type: str,
        method: str = "pooled_t",
    ):
        """两独立样本连续变量。

        差异性：method 为 pooled_t（合并方差 t / 正态近似迭代）、welch（Welch t）、wilcoxon_mw（Mann-Whitney，ARE≈π/3）。
        优效/非劣/等效：两样本 t + 界值或 TOST，固定为合并方差 t 分位数迭代（忽略 welch/wilcoxon）。
        """
        if sd_t <= 0 or sd_c <= 0:
            raise ValueError("标准差必须大于0")

        r = allocation_ratio
        if r <= 0:
            raise ValueError("分配比例须为正。")

        if test_type == "差异性检验":
            effect = abs(mean_t - mean_c)
        elif test_type == "优效性检验":
            effect = (mean_t - mean_c) - margin
        elif test_type == "非劣效性检验":
            effect = (mean_t - mean_c) + margin
        elif test_type == "等效性检验":
            effect = margin - abs(mean_t - mean_c)
        else:
            raise ValueError(f"不支持的检验类型: {test_type}")

        if effect <= 0:
            raise ValueError(f"预期的均值差异无法满足【{test_type}】的统计学前提。")

        equiv = test_type == "等效性检验"
        m = (method or "pooled_t").lower()
        if m in ("z", "pooled", "pooled_t", "t"):
            m_use = "pooled_t"
        elif m in ("welch", "welch_t"):
            m_use = "welch"
        elif m in ("wilcoxon", "mann_whitney", "wilcoxon_mw", "mw"):
            m_use = "wilcoxon_mw"
        else:
            m_use = "pooled_t"

        if test_type != "差异性检验":
            m_use = "pooled_t"

        def _z_crit() -> tuple[float, float]:
            if test_type == "差异性检验":
                z_a = norm.ppf(1 - alpha / 2)
            elif test_type == "等效性检验":
                z_a = norm.ppf(1 - alpha)
            else:
                z_a = norm.ppf(1 - alpha)
            z_b = norm.ppf(1 - (1 - power) / 2.0) if equiv else norm.ppf(power)
            return z_a, z_b

        z_a, z_b = _z_crit()
        variance_term = (sd_t**2) + ((sd_c**2) / r)
        nt_z = ((z_a + z_b) ** 2) * variance_term / (effect**2)
        nt = max(2, int(math.ceil(nt_z)))

        def _se_pooled(nt_i: int, nc_i: int) -> tuple[float, float]:
            df_p = nt_i + nc_i - 2
            if df_p <= 0:
                sp2 = (sd_t**2 + sd_c**2) / 2.0
            else:
                sp2 = ((nt_i - 1) * sd_t**2 + (nc_i - 1) * sd_c**2) / df_p
            se = math.sqrt(sp2 * (1.0 / nt_i + 1.0 / nc_i))
            return se, float(df_p)

        def _se_welch(nt_i: int, nc_i: int) -> tuple[float, float]:
            se = math.sqrt(sd_t**2 / nt_i + sd_c**2 / nc_i)
            dfw = _welch_df(sd_t**2, sd_c**2, float(nt_i), float(nc_i))
            return se, dfw

        for _ in range(250):
            nc = max(2, int(math.ceil(nt * r)))
            if m_use == "welch":
                se, df = _se_welch(nt, nc)
            else:
                se, df = _se_pooled(nt, nc)
       
            if test_type == "差异性检验":
                t_a = t_dist.ppf(1 - alpha / 2, df)
            elif test_type == "等效性检验":
                t_a = t_dist.ppf(1 - alpha, df)
            else:
                t_a = t_dist.ppf(1 - alpha, df)
            t_b = (
                t_dist.ppf(1 - (1 - power) / 2.0, df)
                if equiv
                else t_dist.ppf(power, df)
            )

            nt_new = max(2.0, ((t_a + t_b) ** 2) * (se**2) / (effect**2))
            nt_next = int(math.ceil(nt_new))
            if nt_next == nt:
                break
            nt = nt_next
        else:
            nt = max(2, int(math.ceil(nt_z)))

        nc = max(2, int(math.ceil(nt * r)))
        nt_ceil, nc_ceil = nt, nc

        if m_use == "wilcoxon_mw":
            scale = math.pi / 3.0
            nt_ceil = max(2, int(math.ceil(nt * scale)))
            nc_ceil = max(2, int(math.ceil(nc * scale)))

        if m_use == "welch":
            method_note = "两独立样本：Welch t 检验（方差不齐；分位数迭代）"
        elif m_use == "wilcoxon_mw":
            method_note = (
                "两独立样本：Mann-Whitney / Wilcoxon 秩和（在合并方差 t 样本量上乘 ARE≈π/3）"
            )
        elif test_type == "等效性检验":
            method_note = "两独立样本：等效性 TOST（两个单侧 t，分位数迭代）"
        elif test_type in ("优效性检验", "非劣效性检验"):
            method_note = "两独立样本：均值差 t 检验 + 界值 Δ（单侧；分位数迭代）"
        else:
            method_note = "两独立样本：合并方差 t 检验（正态近似；分位数迭代）"

        nt_drop = int(math.ceil(nt_ceil / (1 - dropout_rate)))
        nc_drop = int(math.ceil(nc_ceil / (1 - dropout_rate)))

        return {
            "n_treatment": nt_ceil,
            "n_control": nc_ceil,
            "total_sample_size": nt_ceil + nc_ceil,
            "n_treatment_with_dropout": nt_drop,
            "n_control_with_dropout": nc_drop,
            "total_sample_size_with_dropout": nt_drop + nc_drop,
            "method": method_note,
            "two_arm_method": m_use,
        }

    def compute_single_arm_precision(
        self,
        sd: float,
        half_width: float,
        alpha: float,
        dropout_rate: float,
    ) -> dict:
        """单臂均值精度：双侧 CI 半宽 ≤ w，采用 t_{α/2, n−1}·SD/√n（t 分布半宽法，非 z）。"""
        if sd <= 0 or half_width <= 0:
            raise ValueError("标准差与均值 CI 半宽须为正。")
        n = max(2, int(math.ceil((norm.ppf(1 - alpha / 2) * sd / half_width) ** 2)))
        for _ in range(500000):
            df = n - 1
            t_crit = t_dist.ppf(1 - alpha / 2, df)
            margin = t_crit * sd / math.sqrt(n)
            if margin <= half_width:
                n_ceil = n
                break
            n += 1
        else:
            raise ValueError("在合理范围内未能找到满足精度要求的样本量。")

        n_drop = int(math.ceil(n_ceil / (1 - dropout_rate)))
        return {
            "n_treatment": n_ceil,
            "n_control": 0,
            "total_sample_size": n_ceil,
            "n_treatment_with_dropout": n_drop,
            "n_control_with_dropout": 0,
            "total_sample_size_with_dropout": n_drop,
            "method": "单臂均值精度：t 分布半宽法（双侧 CI）",
            "precision": {"half_width": half_width, "alpha": alpha},
        }

    def compute_single_arm(
        self,
        mean_t: float,
        mean_ref: float,
        sd: float,
        alpha: float,
        power: float,
        margin: float,
        dropout_rate: float,
        test_type: str,
        method: str = "t",
    ):
        """单臂连续变量假设检验：优效/非劣/等效均用单样本 t（分位数迭代）；差异性可选 t 或 Wilcoxon。

        method: 仅「差异性检验」有效 —— \"t\"（单样本 t）或 \"wilcoxon\"（符号秩，ARE 校正）；其余检验类型固定为 t。
        """
        if sd <= 0:
            raise ValueError("标准差必须大于0")
        m = (method or "t").lower()
        if m not in ("t", "wilcoxon", "z"):
            raise ValueError("单臂差异性方法须为 t 或 wilcoxon（z 已弃用，将按 t 处理）")
        if m == "z":
            m = "t"

        if test_type == "差异性检验":
            effect = abs(mean_t - mean_ref)
        elif test_type == "优效性检验":
            effect = (mean_t - mean_ref) - margin
        elif test_type == "非劣效性检验":
            effect = (mean_t - mean_ref) + margin
        elif test_type == "等效性检验":
            effect = margin - abs(mean_t - mean_ref)
        else:
            raise ValueError(f"不支持的检验类型: {test_type}")

        if effect <= 0:
            raise ValueError(f"预期的均值差异无法满足【{test_type}】的统计学前提。")

        equiv_power = test_type == "等效性检验"

        def _n_z() -> float:
            if test_type == "差异性检验":
                z_a = norm.ppf(1 - alpha / 2)
            elif test_type == "等效性检验":
                z_a = norm.ppf(1 - alpha)
            else:
                z_a = norm.ppf(1 - alpha)
            z_b = norm.ppf(1 - (1 - power) / 2.0) if equiv_power else norm.ppf(power)
            return ((z_a + z_b) ** 2) * (sd ** 2) / (effect ** 2)

        def _n_t_iter() -> int:
            """单样本 t：t 分位数迭代。"""
            n = max(2, int(math.ceil(_n_z())))
            for _ in range(200):
                df = n - 1
                if test_type == "差异性检验":
                    t_a = t_dist.ppf(1 - alpha / 2, df)
                elif test_type == "等效性检验":
                    t_a = t_dist.ppf(1 - alpha, df)
                else:
                    t_a = t_dist.ppf(1 - alpha, df)
                t_b = (
                    t_dist.ppf(1 - (1 - power) / 2.0, df)
                    if equiv_power
                    else t_dist.ppf(power, df)
                )
                n_new = max(2.0, ((t_a + t_b) ** 2) * (sd ** 2) / (effect ** 2))
                n_next = int(math.ceil(n_new))
                if n_next == n:
                    return n
                n = n_next
            return n

        n_t = _n_t_iter()
        if test_type != "差异性检验":
            n_ceil = n_t
            method_note = "单臂：单样本 t 检验（分位数迭代；优效/非劣/等效）"
            m_out = "t"
        elif m == "wilcoxon":
            n_ceil = int(math.ceil(n_t * (math.pi / 3)))
            method_note = (
                "单臂：Wilcoxon 符号秩（在单样本 t 样本量上乘 ARE=π/3≈1.047）"
            )
            m_out = "wilcoxon"
        else:
            n_ceil = n_t
            method_note = "单臂：单样本 t 检验（分位数迭代；差异性）"
            m_out = "t"

        n_drop = int(math.ceil(n_ceil / (1 - dropout_rate)))
        return {
            "n_treatment": n_ceil,
            "n_control": 0,
            "total_sample_size": n_ceil,
            "n_treatment_with_dropout": n_drop,
            "n_control_with_dropout": 0,
            "total_sample_size_with_dropout": n_drop,
            "method": method_note,
            "single_arm_method": m_out,
        }

    def compute_paired_precision(
        self,
        sd_diff: float,
        half_width: float,
        alpha: float,
        dropout_rate: float,
    ) -> dict:
        """配对差值均值精度：双侧 CI 半宽 ≤ w，采用 t_{α/2, n−1}·SD_diff/√n（与单臂同形，自由度为对子数减 1）。"""
        if sd_diff <= 0 or half_width <= 0:
            raise ValueError("配对差值标准差与 CI 半宽须为正。")
        n = max(2, int(math.ceil((norm.ppf(1 - alpha / 2) * sd_diff / half_width) ** 2)))
        for _ in range(500000):
            df = n - 1
            t_crit = t_dist.ppf(1 - alpha / 2, df)
            margin = t_crit * sd_diff / math.sqrt(n)
            if margin <= half_width:
                n_ceil = n
                break
            n += 1
        else:
            raise ValueError("在合理范围内未能找到满足精度要求的样本量。")

        n_drop = int(math.ceil(n_ceil / (1 - dropout_rate)))
        return {
            "n_treatment": n_ceil,
            "n_control": 0,
            "total_sample_size": n_ceil,
            "n_treatment_with_dropout": n_drop,
            "n_control_with_dropout": 0,
            "total_sample_size_with_dropout": n_drop,
            "method": "配对差值均值精度：t 分布半宽法（双侧 CI，对子数 n）",
            "precision": {"half_width": half_width, "alpha": alpha},
        }

    def compute_paired(
        self,
        mean_diff: float,
        sd_diff: float,
        alpha: float,
        power: float,
        margin: float,
        dropout_rate: float,
        test_type: str,
        method: str = "t",
    ):
        """配对两组连续变量：基于配对差值；配对 t 分位数迭代（非单纯 z）。

        method: 仅「差异性检验」有效 —— \"t\"（配对 t）或 \"wilcoxon\"（符号秩 ARE≈π/3）；优效/非劣/等效固定配对 t（等效用 TOST 思路）。
        """
        if sd_diff <= 0:
            raise ValueError("配对差值标准差必须大于0")

        m = (method or "t").lower()
        if m not in ("t", "wilcoxon", "z"):
            raise ValueError("配对差异性方法须为 t 或 wilcoxon（z 已弃用，将按 t 处理）")
        if m == "z":
            m = "t"

        if test_type == "差异性检验":
            effect = abs(mean_diff)
        elif test_type == "优效性检验":
            effect = mean_diff - margin
        elif test_type == "非劣效性检验":
            effect = mean_diff + margin
        elif test_type == "等效性检验":
            effect = margin - abs(mean_diff)
        else:
            raise ValueError(f"不支持的检验类型: {test_type}")

        if effect <= 0:
            raise ValueError(f"预期的配对差异无法满足【{test_type}】的统计学前提。")

        equiv_power = test_type == "等效性检验"

        def _n_z() -> float:
            if test_type == "差异性检验":
                z_a = norm.ppf(1 - alpha / 2)
            elif test_type == "等效性检验":
                z_a = norm.ppf(1 - alpha)
            else:
                z_a = norm.ppf(1 - alpha)
            z_b = norm.ppf(1 - (1 - power) / 2.0) if equiv_power else norm.ppf(power)
            return ((z_a + z_b) ** 2) * (sd_diff ** 2) / (effect ** 2)

        def _n_t_iter() -> int:
            n = max(2, int(math.ceil(_n_z())))
            for _ in range(200):
                df = n - 1
                if test_type == "差异性检验":
                    t_a = t_dist.ppf(1 - alpha / 2, df)
                elif test_type == "等效性检验":
                    t_a = t_dist.ppf(1 - alpha, df)
                else:
                    t_a = t_dist.ppf(1 - alpha, df)
                t_b = (
                    t_dist.ppf(1 - (1 - power) / 2.0, df)
                    if equiv_power
                    else t_dist.ppf(power, df)
                )
                n_new = max(2.0, ((t_a + t_b) ** 2) * (sd_diff ** 2) / (effect ** 2))
                n_next = int(math.ceil(n_new))
                if n_next == n:
                    return n
                n = n_next
            return n

        n_t = _n_t_iter()
        if test_type != "差异性检验":
            n_ceil = n_t
            method_note = "配对两组：配对 t 检验（分位数迭代；优效/非劣/等效）"
            m_out = "t"
        elif m == "wilcoxon":
            n_ceil = int(math.ceil(n_t * (math.pi / 3)))
            method_note = (
                "配对两组：Wilcoxon 符号秩（在配对 t 样本量上乘 ARE=π/3≈1.047）"
            )
            m_out = "wilcoxon"
        else:
            n_ceil = n_t
            method_note = "配对两组：配对 t 检验（分位数迭代；差异性）"
            m_out = "t"

        if test_type == "等效性检验":
            method_note = "配对两组：等效性 TOST（两个单侧配对 t，分位数迭代）"
        elif test_type in ("优效性检验", "非劣效性检验"):
            method_note = "配对两组：配对均值差 t 检验 + 界值 Δ（单侧；分位数迭代）"

        n_drop = int(math.ceil(n_ceil / (1 - dropout_rate)))

        return {
            "n_treatment": n_ceil,
            "n_control": 0,
            "total_sample_size": n_ceil,
            "n_treatment_with_dropout": n_drop,
            "n_control_with_dropout": 0,
            "total_sample_size_with_dropout": n_drop,
            "method": method_note,
            "paired_method": m_out,
        }

    def compute_three_arm_anova(
        self,
        mean1: float,
        mean2: float,
        mean3: float,
        sd_within: float,
        alpha: float,
        power: float,
        dropout_rate: float,
    ) -> dict:
        """三组独立样本均值差异性检验样本量（单因素方差分析 F 检验，等额分配）

        基于非中心 F 分布；组内方差齐性假设下 sd_within 为共同标准差。
        """
        if sd_within <= 0:
            raise ValueError("组内标准差须大于 0。")
        k = 3
        grand = (mean1 + mean2 + mean3) / 3.0
        var_between = (
            (mean1 - grand) ** 2 + (mean2 - grand) ** 2 + (mean3 - grand) ** 2
        ) / k
        if var_between <= 1e-12:
            raise ValueError("三组均值相同，无法计算样本量。")
        f_effect_sq = var_between / (sd_within ** 2)
        df1 = k - 1
        n = 5
        while n < 500000:
            df2 = k * n - k
            lambda_nc = n * k * f_effect_sq
            f_crit = f_dist.ppf(1 - alpha, df1, df2)
            pwr = 1 - ncf.cdf(f_crit, df1, df2, lambda_nc)
            if pwr >= power:
                break
            n += 1
        else:
            raise ValueError("在合理范围内未能找到满足检验效能的样本量，请检查参数。")

        n_per_group_drop = int(math.ceil(n / (1 - dropout_rate)))
        return {
            "n_per_group": n,
            "total_sample_size": k * n,
            "cohens_f": math.sqrt(f_effect_sq),
            "n_per_group_with_dropout": n_per_group_drop,
            "total_sample_size_with_dropout": k * n_per_group_drop,
            "method": "三组独立样本：单因素 ANOVA（F 检验；非中心 F，等额分配；方差齐性）",
            "three_arm_method": "anova_f",
        }

    def _two_sample_equal_n_pooled(
        self,
        delta: float,
        sd: float,
        alpha: float,
        power: float,
    ) -> int:
        """两独立样本、等额分配、合并方差 t 检验所需每组样本量（双侧）。"""
        if delta <= 0 or sd <= 0:
            raise ValueError("效应量与标准差须为正。")
        z_a = norm.ppf(1 - alpha / 2)
        z_b = norm.ppf(power)
        n = max(2, int(math.ceil(2 * (z_a + z_b) ** 2 * sd**2 / delta**2)))
        for _ in range(250):
            df = 2 * n - 2
            t_a = t_dist.ppf(1 - alpha / 2, df)
            t_b = t_dist.ppf(power, df)
            n_new = max(2.0, 2 * (t_a + t_b) ** 2 * sd**2 / delta**2)
            n_next = int(math.ceil(n_new))
            if n_next == n:
                return n
            n = n_next
        return n

    def compute_three_arm_kruskal_wallis(
        self,
        mean1: float,
        mean2: float,
        mean3: float,
        sd_within: float,
        alpha: float,
        power: float,
        dropout_rate: float,
    ) -> dict:
        """在单因素 ANOVA 样本量基础上，按正态假定下 K-W 相对 F 的渐近相对效率 ARE≈0.955 放大每组 n。"""
        base = self.compute_three_arm_anova(
            mean1, mean2, mean3, sd_within, alpha, power, dropout_rate
        )
        n_ano = base["n_per_group"]
        n_kw = int(math.ceil(n_ano / 0.955))
        n_per_group_drop = int(math.ceil(n_kw / (1 - dropout_rate)))
        return {
            "n_per_group": n_kw,
            "total_sample_size": 3 * n_kw,
            "cohens_f": base.get("cohens_f"),
            "n_per_group_with_dropout": n_per_group_drop,
            "total_sample_size_with_dropout": 3 * n_per_group_drop,
            "method": (
                "三组独立样本：Kruskal-Wallis（在 ANOVA 每组样本量上按 ARE≈0.955 放大；"
                "正态假定下近似）"
            ),
            "three_arm_method": "kw",
            "anova_n_per_group_reference": n_ano,
        }

    def compute_three_arm_linear_trend(
        self,
        mean1: float,
        mean2: float,
        mean3: float,
        sd_within: float,
        alpha: float,
        power: float,
        dropout_rate: float,
    ) -> dict:
        """三组等距水平下的线性趋势（正交对比 −1,0,+1）：双侧检验，非中心 t 近似（误差 df=3n−3）。"""
        if sd_within <= 0:
            raise ValueError("组内标准差须大于 0。")
        delta = abs(mean3 - mean1)
        if delta <= 1e-12:
            raise ValueError(
                "有序线性趋势（对比 −1,0,+1）要求组1与组3预期均值不同；请调整均值或换用 ANOVA。"
            )
        _ = mean2  # 参与方案描述；对比检验量主要反映 μ3−μ1

        n = 2
        while n < 500000:
            df = 3 * n - 3
            ncp = (delta / sd_within) * math.sqrt(n / 2.0)
            t_crit = t_dist.ppf(1 - alpha / 2, df)
            between = nct_dist.cdf(t_crit, df, ncp) - nct_dist.cdf(-t_crit, df, ncp)
            pwr = 1.0 - between
            if pwr >= power:
                break
            n += 1
        else:
            raise ValueError("在合理范围内未能找到满足检验效能的样本量，请检查参数。")

        n_per_group_drop = int(math.ceil(n / (1 - dropout_rate)))
        return {
            "n_per_group": n,
            "total_sample_size": 3 * n,
            "cohens_f": None,
            "n_per_group_with_dropout": n_per_group_drop,
            "total_sample_size_with_dropout": 3 * n_per_group_drop,
            "method": (
                "三组独立样本：有序线性趋势（ANOVA 线性对比 −1,0,+1；非中心 t，双侧；"
                "Jonckheere–Terpstra 等非参趋势检验样本量常与此接近，可按 ARE 微调）"
            ),
            "three_arm_method": "linear_trend",
        }

    def compute_three_arm_dunnett(
        self,
        mean1: float,
        mean2: float,
        mean3: float,
        control_index: int,
        sd_within: float,
        alpha: float,
        power: float,
        dropout_rate: float,
    ) -> dict:
        """两组试验臂分别与同一对照比较（Dunnett）；样本量按 Bonferroni（每对比 α/2，双侧）保守近似两两合并方差 t。"""
        if sd_within <= 0:
            raise ValueError("组内标准差须大于 0。")
        if control_index not in (0, 1, 2):
            raise ValueError("对照组须为 0、1 或 2（对应组1–组3）。")
        means = [mean1, mean2, mean3]
        mc = means[control_index]
        trt = [means[i] for i in range(3) if i != control_index]
        deltas = [abs(t - mc) for t in trt]
        delta = min(deltas)
        if delta <= 1e-12:
            raise ValueError(
                "Dunnett 需至少一个试验组与对照的预期均值差非零；请调整均值或对照组选择。"
            )
        # 两个试验组 vs 对照：FWER α 下 Bonferroni 每对比 α/2（双侧）
        alpha_adj = alpha / 2.0
        n = self._two_sample_equal_n_pooled(delta, sd_within, alpha_adj, power)
        n_per_group_drop = int(math.ceil(n / (1 - dropout_rate)))
        return {
            "n_per_group": n,
            "total_sample_size": 3 * n,
            "cohens_f": None,
            "n_per_group_with_dropout": n_per_group_drop,
            "total_sample_size_with_dropout": 3 * n_per_group_drop,
            "method": (
                "三组独立样本：Dunnett（多试验组 vs 对照；Bonferroni α/2 保守近似合并方差 t，"
                "等额分配；正式分析可用 Dunnett 精确界值复核）"
            ),
            "three_arm_method": "dunnett",
            "dunnett_min_pairwise_delta": delta,
        }