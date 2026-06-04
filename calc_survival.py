# calc_survival.py
"""生存终点样本量：Log-rank（Schoenfeld / Freedman）、Lachin–Foulkes 事件概率、里程碑生存率。"""

from __future__ import annotations

import math
from typing import Optional

from scipy.stats import norm


def _hr_from_medians(median_t: float, median_c: float) -> float:
    """指数分布下 HR = h_T/h_C = λ_T/λ_C = m_C/m_T（试验组为 T，对照为 C）。"""
    if median_t <= 0 or median_c <= 0:
        raise ValueError("中位生存期须大于 0。")
    return median_c / median_t


def _lf_event_prob_exponential(median: float, T_a: float, T_f: float) -> float:
    """均匀入组 [0,T_a]、最短随访 T_f、指数生存；单臂在观察结束前的累积事件概率（Lachin–Foulkes 型）。"""
    if median <= 0:
        raise ValueError("中位生存期须大于 0。")
    if T_a <= 0:
        raise ValueError("入组期须大于 0。")
    if T_f < 0:
        raise ValueError("最短随访期须 ≥ 0。")
    lam = math.log(2.0) / median
    # ∫_{T_f}^{T_a+T_f} (1 - e^{-λ s}) ds / T_a
    ta_tf = T_a + T_f
    integ = (ta_tf - T_f) + (1.0 / lam) * (
        math.exp(-lam * ta_tf) - math.exp(-lam * T_f)
    )
    return integ / T_a


class SurvivalCalculator:
    def resolve_hr_two_arm(
        self,
        hr_input_mode: str,
        hr_direct: Optional[float],
        median_t: Optional[float],
        median_c: Optional[float],
    ) -> float:
        if hr_input_mode == "direct":
            if hr_direct is None or hr_direct <= 0:
                raise ValueError("直接输入的 HR 须为正。")
            return float(hr_direct)
        if median_t is None or median_c is None:
            raise ValueError("中位数推算 HR 时须填写试验组与对照组中位生存期。")
        return _hr_from_medians(median_t, median_c)

    def resolve_hr_single_arm(
        self,
        hr_input_mode: str,
        hr_direct: Optional[float],
        median_t: Optional[float],
        median_ref: Optional[float],
    ) -> float:
        if hr_input_mode == "direct":
            if hr_direct is None or hr_direct <= 0:
                raise ValueError("直接输入的 HR 须为正。")
            return float(hr_direct)
        if median_t is None or median_ref is None:
            raise ValueError("中位数推算 HR 时须填写试验组与参考中位生存期。")
        return _hr_from_medians(median_t, median_ref)

    def _z_for_test(self, test_type: str, alpha: float, power: float) -> tuple[float, float]:
        """双侧/单侧 α 与 β；等效性用 TOST 近似 (2 z_α + z_β)² 形式。"""
        if test_type == "差异性检验":
            z_a = norm.ppf(1 - alpha / 2)
            z_b = norm.ppf(power)
            return z_a, z_b
        if test_type in ("优效性检验", "非劣效性检验"):
            z_a = norm.ppf(1 - alpha)
            z_b = norm.ppf(power)
            return z_a, z_b
        if test_type == "等效性检验":
            z_a = norm.ppf(1 - alpha)
            z_b = norm.ppf(1 - (1 - power) / 2.0)
            return z_a, z_b
        raise ValueError(f"不支持的生存检验类型: {test_type}")

    def _log_hr_effect_two_arm(
        self,
        test_type: str,
        hr: float,
        margin_hr: Optional[float],
        equiv_delta: Optional[float],
    ) -> float:
        """返回 |ln HR| 或等价于 log-rank 的效应（用于分母）。"""
        if test_type == "差异性检验":
            if abs(math.log(hr)) < 1e-12:
                raise ValueError("预期 HR 过于接近 1，差异性检验无法估计样本量。")
            return abs(math.log(hr))
        if test_type == "优效性检验":
            if hr >= 1.0:
                raise ValueError("优效性（试验更优）要求预期 HR < 1。")
            return abs(math.log(hr))
        if test_type == "非劣效性检验":
            if margin_hr is None or margin_hr <= 1.0:
                raise ValueError("非劣效须指定非劣界值 HR₀（>1，如 1.3）。")
            if hr >= margin_hr:
                raise ValueError("预期 HR 须小于非劣界值 HR₀。")
            return math.log(margin_hr) - math.log(hr)
        if test_type == "等效性检验":
            if equiv_delta is None or equiv_delta <= 1.0:
                raise ValueError("等效性须指定 HR 等效界值 Δ（>1，如 1.25）。")
            if hr <= 0 or hr < 1.0 / equiv_delta or hr > equiv_delta:
                raise ValueError("预期 HR 须落在等效区间 (1/Δ, Δ) 内。")
            return math.log(equiv_delta)
        raise ValueError(f"不支持的检验类型: {test_type}")

    def compute_two_arm(
        self,
        *,
        hr_input_mode: str,
        hr_direct: Optional[float],
        median_t: Optional[float],
        median_c: Optional[float],
        alpha: float,
        power: float,
        allocation_ratio: float,
        dropout_rate: float,
        test_type: str,
        margin_hr: Optional[float],
        equiv_delta: Optional[float],
        event_mode: str,
        event_rate: Optional[float],
        T_a: Optional[float],
        T_f: Optional[float],
        lr_method: str,
    ) -> dict:
        """双臂：Log-rank 所需总事件数 → 总样本量；事件概率可直接填或用 Lachin–Foulkes（指数、两组分别算）。"""
        hr = self.resolve_hr_two_arm(hr_input_mode, hr_direct, median_t, median_c)
        r = allocation_ratio
        if r <= 0:
            raise ValueError("分配比例须为正。")
        p_t = 1.0 / (1.0 + r)
        p_c = r / (1.0 + r)

        z_a, z_b = self._z_for_test(test_type, alpha, power)
        if test_type == "等效性检验":
            z_sum = 2 * z_a + z_b
            log_eff = self._log_hr_effect_two_arm(test_type, hr, margin_hr, equiv_delta)
            events = int(math.ceil((z_sum**2) / (p_t * p_c * (log_eff**2))))
        else:
            log_eff = self._log_hr_effect_two_arm(test_type, hr, margin_hr, equiv_delta)
            if log_eff <= 0:
                raise ValueError("效应量（对数 HR）须为正。")
            if lr_method == "freedman":
                e_s = ((z_a + z_b) ** 2) / (p_t * p_c * (log_eff**2))
                fac = (hr + 1.0) ** 2 / (4.0 * hr)
                events = int(math.ceil(e_s * fac))
            else:
                events = int(
                    math.ceil(((z_a + z_b) ** 2) / (p_t * p_c * (log_eff**2)))
                )

        if event_mode == "direct":
            if event_rate is None or not (0 < event_rate <= 1):
                raise ValueError("直接模式须填写 (0,1] 内的总体预期事件概率。")
            d_bar = event_rate
        else:
            if median_t is None or median_c is None:
                raise ValueError("Lachin–Foulkes 模式须填写两组中位生存期以推算各组事件概率。")
            if T_a is None or T_f is None or T_a <= 0 or T_f < 0:
                raise ValueError("须填写入组期与最短随访期（月）。")
            d_t = _lf_event_prob_exponential(median_t, T_a, T_f)
            d_c = _lf_event_prob_exponential(median_c, T_a, T_f)
            d_bar = p_t * d_t + p_c * d_c
            if d_bar <= 0:
                raise ValueError("推算的总体事件概率须大于 0。")

        total_n = int(math.ceil(events / d_bar))
        nt_ceil = int(math.ceil(total_n * p_t))
        nc_ceil = int(math.ceil(total_n * p_c))
        nt_drop = int(math.ceil(nt_ceil / (1 - dropout_rate)))
        nc_drop = int(math.ceil(nc_ceil / (1 - dropout_rate)))

        method = (
            "双臂：Log-rank（Schoenfeld）"
            if lr_method == "schoenfeld"
            else "双臂：Log-rank（Freedman 校正）"
        )
        if test_type == "等效性检验":
            method = "双臂：Log-rank（Schoenfeld）；等效性 TOST 近似"
        if event_mode == "lachin":
            method += "；事件概率 Lachin–Foulkes（均匀入组+指数生存）"

        return {
            "hr": hr,
            "events": events,
            "n_treatment": nt_ceil,
            "n_control": nc_ceil,
            "total_sample_size": nt_ceil + nc_ceil,
            "n_treatment_with_dropout": nt_drop,
            "n_control_with_dropout": nc_drop,
            "total_sample_size_with_dropout": nt_drop + nc_drop,
            "method": method,
            "event_prob_overall": d_bar,
        }

    def _log_hr_effect_single_arm(
        self,
        test_type: str,
        hr: float,
        margin_hr: Optional[float],
        equiv_delta: Optional[float],
    ) -> float:
        if test_type == "差异性检验":
            if abs(math.log(hr)) < 1e-12:
                raise ValueError("预期 HR 过于接近 1。")
            return abs(math.log(hr))
        if test_type == "优效性检验":
            if hr >= 1.0:
                raise ValueError("优效性要求预期 HR < 1（相对历史对照）。")
            return abs(math.log(hr))
        if test_type == "非劣效性检验":
            if margin_hr is None or margin_hr <= 1.0:
                raise ValueError("非劣效须指定非劣界值 HR₀（>1）。")
            if hr >= margin_hr:
                raise ValueError("预期 HR 须小于 HR₀。")
            return math.log(margin_hr) - math.log(hr)
        if test_type == "等效性检验":
            if equiv_delta is None or equiv_delta <= 1.0:
                raise ValueError("等效性须指定界值 Δ（>1）。")
            if hr < 1.0 / equiv_delta or hr > equiv_delta:
                raise ValueError("预期 HR 须落在 (1/Δ, Δ) 内。")
            return math.log(equiv_delta)
        raise ValueError(f"不支持的检验类型: {test_type}")

    def compute_single_arm(
        self,
        *,
        hr_input_mode: str,
        hr_direct: Optional[float],
        median_t: Optional[float],
        median_ref: Optional[float],
        alpha: float,
        power: float,
        dropout_rate: float,
        test_type: str,
        margin_hr: Optional[float],
        equiv_delta: Optional[float],
        event_mode: str,
        event_rate: Optional[float],
        T_a: Optional[float],
        T_f: Optional[float],
    ) -> dict:
        """单臂：相对历史对照的 HR（Log-rank 事件数近似）。"""
        hr = self.resolve_hr_single_arm(hr_input_mode, hr_direct, median_t, median_ref)
        z_a, z_b = self._z_for_test(test_type, alpha, power)
        if test_type == "等效性检验":
            z_sum = 2 * z_a + z_b
            log_eff = self._log_hr_effect_single_arm(test_type, hr, margin_hr, equiv_delta)
            events = int(math.ceil((z_sum**2) / (log_eff**2)))
        else:
            log_eff = self._log_hr_effect_single_arm(test_type, hr, margin_hr, equiv_delta)
            if log_eff <= 0:
                raise ValueError("效应量须为正。")
            events = int(math.ceil(((z_a + z_b) ** 2) / (log_eff**2)))

        if event_mode == "direct":
            if event_rate is None or not (0 < event_rate <= 1):
                raise ValueError("须填写 (0,1] 内的预期事件概率。")
            d_bar = event_rate
        else:
            if median_t is None:
                raise ValueError("Lachin–Foulkes 须填写试验组中位生存期（月）。")
            if T_a is None or T_f is None:
                raise ValueError("须填写入组期与最短随访期。")
            d_bar = _lf_event_prob_exponential(median_t, T_a, T_f)

        total_n = int(math.ceil(events / d_bar))
        total_n_drop = int(math.ceil(total_n / (1 - dropout_rate)))
        method = "单臂：相对历史对照的 Log-rank 近似（事件数）"
        if test_type == "等效性检验":
            method += "；等效 TOST 近似"
        if event_mode == "lachin":
            method += "；事件概率 Lachin–Foulkes"
        return {
            "hr": hr,
            "events": events,
            "n_treatment": total_n,
            "n_control": 0,
            "total_sample_size": total_n,
            "n_treatment_with_dropout": total_n_drop,
            "n_control_with_dropout": 0,
            "total_sample_size_with_dropout": total_n_drop,
            "method": method,
            "event_prob_overall": d_bar,
        }

    def compute_single_arm_milestone(
        self,
        *,
        t_months: float,
        s0: float,
        s1: float,
        alpha: float,
        power: float,
        dropout_rate: float,
        test_type: str,
    ) -> dict:
        """固定时间点的生存率：正态近似二样本率差（历史率已知）。"""
        if not (0 < s0 < 1 and 0 < s1 < 1):
            raise ValueError("参考与预期的生存率须在 (0,1) 内。")
        if t_months <= 0:
            raise ValueError("里程碑时间须为正（月）。")
        delta = s1 - s0
        if abs(delta) < 1e-8:
            raise ValueError("预期生存率与参考须不同。")

        if test_type == "差异性检验":
            z_a = norm.ppf(1 - alpha / 2)
            z_b = norm.ppf(power)
        elif test_type in ("优效性检验", "非劣效性检验"):
            z_a = norm.ppf(1 - alpha)
            z_b = norm.ppf(power)
            if test_type == "优效性检验" and s1 <= s0:
                raise ValueError("优效性要求预期生存率高于参考。")
            if test_type == "非劣效性检验" and s1 < s0:
                raise ValueError("非劣效要求预期生存率不低于参考（或请换用 HR 模式设界值）。")
        elif test_type == "等效性检验":
            raise ValueError("里程碑模式暂不支持等效性，请换用 HR 模式或简化假设。")
        else:
            raise ValueError(f"不支持的检验类型: {test_type}")

        # Var(p1 - p0) ≈ p1(1-p1)/n  （单臂估计 p1，p0 已知固定）
        n = ((z_a + z_b) ** 2) * s1 * (1 - s1) / (delta**2)
        n_ceil = int(math.ceil(n))
        n_drop = int(math.ceil(n_ceil / (1 - dropout_rate)))
        return {
            "hr": None,
            "events": None,
            "n_treatment": n_ceil,
            "n_control": 0,
            "total_sample_size": n_ceil,
            "n_treatment_with_dropout": n_drop,
            "n_control_with_dropout": 0,
            "total_sample_size_with_dropout": n_drop,
            "method": f"单臂：里程碑 {t_months:.1f} 月生存率（正态近似；S₀={s0:.3f}，S₁={s1:.3f}）",
            "event_prob_overall": None,
            "milestone_t": t_months,
            "s0": s0,
            "s1": s1,
        }
