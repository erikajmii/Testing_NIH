

//Function that calculates t_core
export function two_nodes_optimized(
    tdb,
    tr,
    v,
    met,
    clo,
    burn_surface,
    length_time_simulation,
    vapor_pressure,
    wme,
    body_surface_area,
    p_atmospheric,
    body_position,
    calculate_ce = false,
    max_skin_blood_flow = 90,
    max_sweating = 500,
    w_max = false
  ) {
    // Initial variables as defined in the ASHRAE 55-2020
    const air_speed = Math.max(v, 0.1);
    const k_clo = 0.25;
    const body_weight = 70; // body weight in kg
    const met_factor = 58.2; // met conversion factor
    const sbc = 0.000000056697; // Stefan-Boltzmann constant (W/m2K4)
    const c_sw = 170; // driving coefficient for regulatory sweating
    const c_dil = 120; // driving coefficient for vasodilation ashrae says 50 see page 9.19
    const c_str = 0.5; // driving coefficient for vasoconstriction
  
    const temp_skin_neutral = 33.7;
    const temp_core_neutral = 36.8;
    let alfa = 0.1;
    const temp_body_neutral = alfa * temp_skin_neutral + (1 - alfa) * temp_core_neutral;
    const skin_blood_flow_neutral = 6.3;
  
    let t_skin = temp_skin_neutral;
    let t_core = temp_core_neutral;
    let m_bl = skin_blood_flow_neutral;

    let e_skin = 0.1 * met;
    let q_sensible = 0;
    let w = 0;
    let _set = 0;
    let e_rsw = 0;
    let e_diff = 0;
    let e_max = 0;
    let m_rsw = 0;
    let q_res = 0;
    let et = 0;
    let e_req = 0;
    let r_ea = 0;
    let r_ecl = 0;
    let c_res = 0;
    let pressure_in_atmospheres = p_atmospheric / 101325;
    let n_simulation = 0;
    let r_clo = 0.155 * clo;
    let f_a_cl = 1.0 + 0.15 * clo;
    let lr = 2.2 / pressure_in_atmospheres;
    let rm = (met - wme) * met_factor;
    let m = met * met_factor;
    let e_comfort = 0.42 * (rm - met_factor);

    if (e_comfort < 0) {
        e_comfort = 0;
    }

    let i_cl = 1.0;

    if (clo > 0) {
        i_cl = 0.45;
    }

    if (!w_max) {
        w_max = 0.38 * Math.pow(air_speed, -0.29);

        if (clo > 0) {
            w_max = 0.59 * Math.pow(air_speed, -0.08);
        }
    }

    w_max *= burn_surface;

    let h_cc = 3.0 * Math.pow(pressure_in_atmospheres, 0.53);
    let h_fc = 8.600001 * Math.pow((air_speed * pressure_in_atmospheres), 0.53);
    h_cc = Math.max(h_cc, h_fc);

    if (!calculate_ce && met > 0.85) {
        let h_c_met = 5.66 * Math.pow((met - 0.85), 0.39);
        h_cc = Math.max(h_cc, h_c_met);
    }

    let h_r = 4.7;
    let h_t = h_r + h_cc;
    let r_a = 1.0 / (f_a_cl * h_t);
    let t_op = (h_r * tr + h_cc * tdb) / h_t;
    let t_body = alfa * t_skin + (1 - alfa) * t_core;
    q_res = 0.0023 * m * (44.0 - vapor_pressure);
    c_res = 0.0014 * m * (34.0 - tdb);

    while (n_simulation < length_time_simulation) {
        n_simulation++;

        let iteration_limit = 150;  //changed
        let t_cl = (r_a * t_skin + r_clo * t_op) / (r_a + r_clo);
        let n_iterations = 0;
        let tc_converged = false;

        while (!tc_converged) {
            let h_r;

            if (body_position === "sitting") {
                h_r = 4.0 * 0.95 * sbc * ((t_cl + tr) / 2.0 + 273.15) ** 3.0 * 0.7;
            } else {
                h_r = 4.0 * 0.95 * sbc * ((t_cl + tr) / 2.0 + 273.15) ** 3.0 * 0.73;
            }

            h_t = h_r + h_cc;
            r_a = 1.0 / (f_a_cl * h_t);
            t_op = (h_r * tr + h_cc * tdb) / h_t;
            let t_cl_new = (r_a * t_skin + r_clo * t_op) / (r_a + r_clo);

            if (Math.abs(t_cl_new - t_cl) <= 0.01) {
                tc_converged = true;
            }

            t_cl = t_cl_new;
            n_iterations++;

            if (n_iterations > iteration_limit) {
                throw new Error("Max iterations exceeded");
            }
        }

        q_sensible = (t_skin - t_op) / (r_a + r_clo);
        let hf_cs = (t_core - t_skin) * (5.28 + 1.163 * m_bl);
        let s_core = m - hf_cs - q_res - c_res - wme;
        let s_skin = hf_cs - q_sensible - e_skin;
        let tc_sk = 0.97 * alfa * body_weight;
        let tc_cr = 0.97 * (1 - alfa) * body_weight;
        let d_t_sk = (s_skin * body_surface_area) / (tc_sk * 60.0);
        let d_t_cr = (s_core * body_surface_area) / (tc_cr * 60.0);
        t_skin = t_skin + d_t_sk;
        t_core = t_core + d_t_cr;
        t_body = alfa * t_skin + (1 - alfa) * t_core;
        let sk_sig = t_skin - temp_skin_neutral;
        let warm_sk = (sk_sig > 0) * sk_sig;
        let colds = ((-1.0 * sk_sig) > 0) * (-1.0 * sk_sig);
        let c_reg_sig = t_core - temp_core_neutral;
        let c_warm = (c_reg_sig > 0) * c_reg_sig;
        let c_cold = ((-1.0 * c_reg_sig) > 0) * (-1.0 * c_reg_sig);
        let bd_sig = t_body - temp_body_neutral;
        let warm_b = (bd_sig > 0) * bd_sig;
        m_bl = (skin_blood_flow_neutral + c_dil * c_warm) / (1 + c_str * colds);

        if (m_bl > max_skin_blood_flow) {
            m_bl = max_skin_blood_flow;
        }

        if (m_bl < 0.5) {
            m_bl = 0.5;
        }

        m_rsw = c_sw * warm_b * Math.exp(warm_sk / 10.7);

        if (m_rsw > max_sweating) {
            m_rsw = max_sweating;
        }

        e_rsw = 0.68 * m_rsw;
        r_ea = 1.0 / (lr * f_a_cl * h_cc);
        r_ecl = r_clo / (lr * i_cl);
        e_req = rm - q_res - c_res - q_sensible;
        e_max = (Math.exp(18.6686 - 4030.183 / (t_skin + 235.0)) - vapor_pressure) / (r_ea + r_ecl);
        let p_rsw = e_rsw / e_max;
        w = 0.06 + 0.94 * p_rsw;
        e_diff = w * e_max - e_rsw;
        
        if (w > w_max) {
            w = w_max;
            p_rsw = w_max / 0.94;
            e_rsw = p_rsw * e_max;
            e_diff = 0.06 * (1.0 - p_rsw) * e_max;
        }
        
        if (e_max < 0) {
            e_diff = 0;
            e_rsw = 0;
            w = w_max;
        }
        
        e_skin = e_rsw + e_diff; // total evaporative heat loss sweating and vapor diffusion
        m_rsw = e_rsw / 0.68; // back calculating the mass of regulatory sweating as a function of e_rsw
        let met_shivering = 19.4 * colds * c_cold; // met shivering W/m2
        m = rm + met_shivering;
        alfa = 0.0417737 + 0.7451833 / (m_bl + 0.585417);
    }

    // Calculate total heat loss from skin
    let q_skin = q_sensible + e_skin;

    // Calculate saturation vapour pressure of water of the skin
    let p_s_sk = Math.exp(18.6686 - 4030.183 / (t_skin + 235.0));

    // Standard environment variables
    let h_r_s = h_r; // standard environment radiative heat transfer coefficient

    let h_c_s = 3.0 * Math.pow(pressure_in_atmospheres, 0.53);
    if (!calculate_ce && met > 0.85) {
        let h_c_met = 5.66 * Math.pow(met - 0.85, 0.39);
        h_c_s = Math.max(h_c_s, h_c_met);
    }
    if (h_c_s < 3.0) {
        h_c_s = 3.0;
    }

    let h_t_s = h_c_s + h_r_s; // sum of convective and radiant heat transfer coefficient W/(m2*K)
    let r_clo_s = 1.52 / ((met - wme / met_factor) + 0.6944) - 0.1835; // thermal resistance of clothing, °C M^2 /W
    let r_cl_s = 0.155 * r_clo_s; // thermal insulation of the clothing in M2K/W
    let f_a_cl_s = 1.0 + k_clo * r_clo_s; // increase in body surface area due to clothing
    let f_cl_s = 1.0 / (1.0 + 0.155 * f_a_cl_s * h_t_s * r_clo_s); // ratio of surface clothed body over nude body
    let i_m_s = 0.45; // permeation efficiency of water vapour through the clothing layer
    let i_cl_s = (i_m_s * h_c_s / h_t_s * (1 - f_cl_s) / (h_c_s / h_t_s - f_cl_s * i_m_s)); // clothing vapor permeation efficiency
    let r_a_s = 1.0 / (f_a_cl_s * h_t_s); // resistance of air layer to dry heat
    let r_ea_s = 1.0 / (lr * f_a_cl_s * h_c_s);
    let r_ecl_s = r_cl_s / (lr * i_cl_s);
    let h_d_s = 1.0 / (r_a_s + r_cl_s);
    let h_e_s = 1.0 / (r_ea_s + r_ecl_s);

    // Calculate Standard Effective Temperature (SET)
    let delta = 0.0001;
    let dx = 100.0;
    let set_old = Math.round(t_skin - q_skin / h_d_s, 2);
    while (Math.abs(dx) > 0.01) {
        let err_1 = q_skin - h_d_s * (t_skin - set_old) - w * h_e_s * (p_s_sk - 0.5 * (Math.exp(18.6686 - 4030.183 / (set_old + 235.0))));
        let err_2 = q_skin - h_d_s * (t_skin - (set_old + delta)) - w * h_e_s * (p_s_sk - 0.5 * (Math.exp(18.6686 - 4030.183 / (set_old + delta + 235.0))));
        let _set = set_old - delta * err_1 / (err_2 - err_1);
        dx = _set - set_old;
        set_old = _set;
    }

    // Calculate Effective Temperature (ET)
    let h_d = 1 / (r_a + r_clo);
    let h_e = 1 / (r_ea + r_ecl);
    let et_old = t_skin - q_skin / h_d;
    delta = 0.0001;
    dx = 100.0;
    while (Math.abs(dx) > 0.01) {
        let err_1 = q_skin - h_d * (t_skin - et_old) - w * h_e * (p_s_sk - 0.5 * (Math.exp(18.6686 - 4030.183 / (et_old + 235.0))));
        let err_2 = q_skin - h_d * (t_skin - (et_old + delta)) - w * h_e * (p_s_sk - 0.5 * (Math.exp(18.6686 - 4030.183 / (et_old + delta + 235.0))));
        let et = et_old - delta * err_1 / (err_2 - err_1);
        dx = et - et_old;
        et_old = et;
    }

    let tbm_l = (0.194 / 58.15) * rm + 36.301; // lower limit for evaporative regulation
    let tbm_h = (0.347 / 58.15) * rm + 36.669; // upper limit for evaporative regulation

    let t_sens = 0.4685 * (t_body - tbm_l); // predicted thermal sensation
    if (t_body >= tbm_l && t_body < tbm_h) {
        t_sens = w_max * 4.7 * (t_body - tbm_l) / (tbm_h - tbm_l);
    } else if (t_body >= tbm_h) {
        t_sens = w_max * 4.7 + 0.4685 * (t_body - tbm_h);
    }

    let disc = 4.7 * (e_rsw - e_comfort) / (e_max * w_max - e_comfort - e_diff); // predicted thermal discomfort
    if (disc <= 0) {
        disc = t_sens;
    }

    // PMV Gagge
    let pmv_gagge = (0.303 * Math.exp(-0.036 * m) + 0.028) * (e_req - e_comfort - e_diff);

    // PMV SET
    let dry_set = h_d_s * (t_skin - _set);
    let e_req_set = rm - c_res - q_res - dry_set;
    let pmv_set = (0.303 * Math.exp(-0.036 * m) + 0.028) * (e_req_set - e_comfort - e_diff);

    // Predicted Percent Satisfied With the Level of Air Movement
    let ps = 100 * (1.13 * Math.sqrt(t_op) - 0.24 * t_op + 2.7 * Math.sqrt(v) - 0.99 * v);

    // Print the thermal sensation and core temperature
    let rounded_t_sens = t_sens.toFixed(2);
    let rounded_t_core = t_core.toFixed(2);

    return rounded_t_core;
}

//Function that casts parameters and calls two_nodes_optimized 
export async function two_nodes(tdb, rh, met, met2, clo, burn_surface, length_time_simulation, wme=0,body_surface_area=1.7,p_atmospheric=101325,body_position="standing",max_skin_blood_flow=90){
    //cast variables
    tdb = parseFloat(tdb);
    let tr = tdb;
    //get windspeed variable but for now keep it
    let v = 0;
    rh = parseFloat(rh);
    met = parseFloat(met2);
    clo = parseFloat(clo);
    burn_surface = (parseFloat(burn_surface));
    burn_surface = (100 - burn_surface) * 0.01; // getting the factor we need to multiply w_max by
    length_time_simulation = parseInt(length_time_simulation);

        const vapor_pressure = (rh * Math.exp(18.6686 - (4030.183 / (tdb + 235.0)))) / 100;

        let t_core = two_nodes_optimized(tdb,tr,v,met,clo,burn_surface,length_time_simulation,vapor_pressure,wme,body_surface_area,p_atmospheric,body_position);
            
        return t_core;
}
