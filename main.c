#include<stdio.h>
#include<math.h>

double Pi = 3.14;
double g = 9.81;
double rho0 = 1.225;

double mtoft = 3.28084;
double mstoknot = 1.94384;
double lbtokg = 0.453592;
double kgtolb = 1 / 0.453592;
double pitorad = 0.0174533;
double rpmtoomega = 0.104719;
double omegatorpm = 1 / 0.104719;

double Get_rho(double rho, double h)
{
	return rho * exp(-0.0296 * h / 304.8);
}

double Get_omega(double rpm)
{
	return rpm * 2 * Pi / 60;
}

double Get_A_wing(double root, double tip, double span)
{
	return (root + tip) * span * 0.5;
}

double Get_vi(double T_rotor, double rho, double A_disk)
{
	return sqrt(T_rotor / (2 * rho * A_disk));
}

double Get_P_i(double T_rotor, double vi)
{
	return T_rotor * vi;
}

double Get_P_0(double sigma, double Cd0, double rho, double A_disk, double v_tip)
{
	return ((sigma * Cd0) / 8) * rho * A_disk * pow(v_tip, 3);
}


int main()
{
	/////////////////
	/// 기?초 정보 ///
	/////////////////

	double W_gross = 3000;
	double W_payload_req = 550;
	double payload_calc = 0.0;

	int i = 0;

	double WL = 168.790948; // kg/m²
	double DL = 100;	// kg/m²
	double AR = 8;
	printf("AR = %lf \n", AR);

	double N_rotor = 4;
	double N_crotor = 3;

	double rpm = 500;

	while (abs(W_payload_req - payload_calc) > 1)
	{
		printf("---------------------------------------\n");
		printf("--------------  %d  ------------------\n", i);
		printf("---------------------------------------\n");
		i++;

		double A_disk = W_gross / (DL * N_rotor);
		printf("A_disk = %lf [m²] \n", A_disk);

		double r_disk = sqrt(A_disk / Pi);
		printf("r_disk = %lf [m] \n", r_disk);

		printf("d_disk = %lf [m]\n", r_disk * 2);

		double Cd0 = 0.04293;

		double eta_m = 0.9;
		double eta_prop = 0.9;

		double V_max = 55.55555556;

		double sigma = 0.2262;

		double Cd_rod = 0.00161;

		double Cd_vt = 0.00099;

		double Cd = 0.00808;

		//double P_hover = sqrt(pow(T_hover, 3) / (2 * rho0 * A_disk * 4));
		//double pHOVER = P_hover / W_gross;
		//printf("P_hover = %lf [W]\n", P_hover);

		//double e = 0.8;
		//double K = pow(Pi * e * AR, -1);
		//double eta_P = 0.8; // 프로펠러 효율.
		//double sigma_air = 11000;
		//double CD0 = 0.01333; // 항공기 전체 항력계수. 계속 수정해야 됨.
		//double pCRUISE = (0.5 * rho0 * pow(V_max, 3) * CD0 * pow(WL, -1) + (2 * K * WL) / (rho0 * sigma_air * V_max)) / (eta_P);

		//double P_cruise = pCRUISE * W_gross;
		//printf("P_cruise = %lf [W]\n", P_cruise);

		//double RoC = 3.333333;
		//double q = 0.5 * rho0 * pow(V_max, 2);
		//double CD_min = 0.001;
		//double pCLIMB = (V_max / eta_P) * (RoC / V_max + q * CD_min / WL + K * WL / q) * 9.80778138;
		//double P_climb = pCLIMB * W_gross;
		//printf("P_climb = %lf [W]\n\n", P_climb);

		//double P_motor = P_hover / 4;

		//double R = 250 * 1000;
		//double eta_m = 0.85;
		//double E = 380;
		//double LD = 6; // 수정 필요.
		//double W_bat = (R * g * W_gross) / (LD * eta_P * eta_m * E * 3600);
		//printf("W_bat = %lf [kg]\n\n", W_bat);

		////////////////////////////////////
		/// Cruise 해석 겸 wing span 찾기 ///
		////////////////////////////////////

		double range_cruise = 50 * 1000; // [m]
		double rho_cruise = Get_rho(rho0, 600);
		double t_cruise = range_cruise / V_max;
		printf("t_cruise = %lf [s]\n", t_cruise);
		double vt_span = 2;
		double vt_root = 1.2;
		double vt_tip = 0.6;
		double A_vt = 0.0;
		double A_wing = W_gross / WL;

		double T1 = 2;
		double alpha = 4 * Pi / 180; // degree to radian
		double L_wing = 0.0;
		double D_wing = 0.0;

		double W_gross_lb = W_gross * kgtolb;


		double fe_fuse = 1.7 * pow(((W_gross_lb) / 1000), 2/3);		// lb
		double D_fuselage = fe_fuse * 0.5 * rho_cruise * pow(V_max, 2);

		double fe_hub = 0.7 * pow(((W_gross_lb) / 1000), 2/3);		// lb
		double D_hub = fe_hub * 0.5 * rho_cruise * pow(V_max, 2);

		double Cl = 0.9313;
		double Cdp = 0.00523;
		double e = 0.8;
		double Cdi = pow(Cl, 2) / (Pi * AR * e);

		double A_rod = Pi * pow(0.1, 2);
		double D_rod = 0.5 * rho_cruise * pow(V_max, 2) * A_rod * Cd_rod;

		A_vt = (vt_root + vt_tip) * vt_span / 2;
		double D_vt = 0.5 * rho_cruise * pow(V_max, 2) * A_vt * Cd_vt;

		while (abs(T1) >= 1)
		{
			L_wing = 0.5 * rho_cruise * pow(V_max, 2) * A_wing * Cl;
			D_wing = 0.5 * rho_cruise * pow(V_max, 2) * A_wing * Cd;

			T1 = W_gross * g - L_wing * cos(alpha) + D_wing * sin(alpha);

			if (abs(T1) < 1)
				break;
			else if ((L_wing * cos(alpha) - D_wing * sin(alpha)) > W_gross * g)
				A_wing = A_wing * 0.5;
			else
				A_wing = A_wing * 1.5;
		}

		WL = W_gross / A_wing;
		//printf("T1 = %lf [N] \n", T1);

		printf("A_wing = %lf [m²]", A_wing);
		printf("AR = %lf \n\n\n", AR);
		//printf("A_wing / 2 = %lf [m²]\n", A_wing / 2);
		printf("WL = %lf [Pa]\n", WL);
		double F_x = D_hub * 4 + D_fuselage + D_rod * 2 + D_vt * 2 + D_wing * cos(alpha) + L_wing * sin(alpha);

		double P_cruise_req = F_x * V_max;
		printf("P_cruise_req = %lf [W] \n", P_cruise_req);

		double P_cruise_prop = (P_cruise_req / N_crotor) / eta_prop;
		printf("P_cruise_prop = %lf [W]\n", P_cruise_prop);

		double P_cruise_motor = P_cruise_prop / eta_m;
		printf("P_cruise_motor = %lf [W]\n", P_cruise_motor);

		double P_cruise_total = P_cruise_motor * N_crotor;
		printf("P_cruise_total = %lf [W]\n", P_cruise_total);

		/////////////
		/// climb ///
		/////////////

		double T = W_gross * g;
		double T_rotor = T / N_rotor;

		double range_climb = 600;
		double omega_climb = Get_omega(rpm);
		double v_tip_climb = omega_climb * r_disk;
		double rho_climb = Get_rho(rho0, 0);
		double v_climb = 3.333333;

		double t_climb = range_climb / v_climb;

		double vi_climb = Get_vi(T_rotor, rho_climb, A_disk);
		double P_i_climb = Get_P_i(T_rotor, vi_climb + v_climb);
		double P_0_climb = Get_P_0(sigma, Cd0, rho_climb, A_disk, v_tip_climb);

		double P_climb_req = (P_i_climb + P_0_climb) * N_rotor;

		printf("P_climb_req = %lf [N]\n", P_climb_req);

		double P_climb_prop = (P_climb_req / N_rotor) / eta_prop;
		printf("P_cilmb_prop = %lf [W]\n", P_climb_prop);

		double P_climb_motor = P_climb_prop / eta_m;
		printf("P_climb_motor = %lf [W]\n", P_climb_motor);

		double P_climb_total = P_climb_motor * N_rotor;
		printf("P_climb_total = %lf [W]\n", P_climb_total);

		/////////////
		/// hover ///
		/////////////

		double t_hover = 0;
		double omega_hover = Get_omega(rpm);
		double v_tip_hover = omega_hover * r_disk;
		double vi = Get_vi(T_rotor, rho_climb, A_disk);

		double P_i_hover = Get_P_i(T_rotor, vi);
		double P_0_hover = Get_P_0(sigma, Cd0, rho0, A_disk, v_tip_hover);

		double P_hover_req = (P_i_hover + P_0_hover) * N_rotor;
		printf("P_hover_req = %lf [W]\n", P_hover_req);

		double P_hover_prop = (P_hover_req / N_rotor) / eta_prop;
		printf("P_hover_prop = %lf [W]\n", P_hover_prop);

		double P_hover_motor = P_hover_prop / eta_m;
		printf("P_hover_motor = %lf [W]\n", P_hover_motor);

		double P_hover_total = P_hover_motor * N_rotor;
		printf("P_hover_total = %lf [W]\n", P_hover_total);

		//while (abs(P_hover_total - P_hover) > 1)
		//{
		//	if (P_hover_total > P_hover)
		//		rpm_hover = rpm_hover * 0.5;
		//	else
		//		rpm_hover = rpm_hover * 1.5;

		//	omega_hover = Get_omega(rpm_hover);
		//	v_tip_hover = omega_hover * r_disk;
		//	P_i_hover = Get_P_i(T_rotor, vi);
		//	P_0_hover = Get_P_0(sigma, Cd0, rho0, A_disk, v_tip_hover);

		//	P_hover_motor = P_i_hover + P_0_hover;
		//	P_hover_total = P_hover_motor * N_rotor;
		//}

		//printf("rpm_hover = %lf\n", rpm_hover);	// 1266

		///////////////
		/// descent ///
		///////////////

		double omega_descent = Get_omega(rpm);
		double v_tip_descent = omega_descent * r_disk;
		double rho_descent = Get_rho(rho0, 600);
		double t_descent = 9 * 60;	// 9 min
		double v_descent = -1.111111;

		double vi_descent = Get_vi(T_rotor, rho_descent, A_disk);
		double P_i_descent = Get_P_i(T_rotor, vi_descent + v_descent);
		double P_0_descent = Get_P_0(sigma, Cd0, rho_descent, A_disk, v_tip_descent);

		double P_descent_req = (P_i_descent + P_0_descent) * N_rotor;
		printf("P_descent_req = %lf [W]\n", P_descent_req);

		double P_descent_prop = (P_descent_req / N_rotor) / eta_prop;
		printf("P_descent_prop = %lf [W]\n", P_descent_prop);

		double P_descent_motor = P_descent_prop / eta_m;
		printf("P_descent_motor = %lf [W]\n", P_descent_motor);

		double E_cruise = P_cruise_total * t_cruise / 3600;
		printf("E_cruise = %lf [Wh]\n\n", E_cruise);

		double E_climb = P_climb_total * (t_climb / 3600);
		printf("E_climb = %lf [Wh]\n\n", E_climb);

		double E_hover = P_hover_total * t_hover / 3600;
		printf("E_hover = %lf [Wh]\n\n", E_hover);

		double E_descent = P_descent_motor * (t_descent / 3600);
		printf("E_descent = %lf [Wh]\n\n", E_descent);

		double E_total = E_hover + E_climb + E_cruise + E_descent;

		double density_bat = 380; // [Wh/kg]
		double W_bat = E_total * 1.2 / density_bat;
		printf("W_bat = %lf [kg]\n", W_bat);

		double P_max_arr[] = { P_cruise_motor, P_climb_motor, P_hover_motor, P_descent_motor };

		int size = sizeof(P_max_arr) / sizeof(P_max_arr[0]);

		// 배열에서 가장 큰 값을 찾기 위한 변수 선언 및 초기화
		double P_max = P_max_arr[0];

		// 배열을 순회하며 가장 큰 값 찾기
		for (int i = 1; i < size; i++)
		{
			if (P_max_arr[i] > P_max)
				P_max = P_max_arr[i];
		}

		printf("P_max = %lf [W]\n", P_max);

		double SP_motor = 10700;
		double W_motor = (N_rotor + N_crotor)*(P_max / SP_motor);
		printf("W_motor = %lf [kg]\n", W_motor);

		//double N = 3.4; // 하중계수.
		//double Rambda = 3.5 * pitorad;
		//double rambda = 1.4 / 1.82519;	// taper ratio. 0.8로 가정
		//double tc = 0.12044;
		//rho_cruise = Get_rho(rho0, 600);
		double V_max_knot = V_max * mstoknot;
		double Ve_knot = 0.5 * rho_cruise * pow(V_max_knot, 2);
		//double A_wing_ft = A_wing * pow(mtoft, 2);

		//double wa = pow(W_gross_lb * N / 100000, 0.65);
		//double wb = pow(AR / cos(Rambda), 0.57);
		//double wc = pow(A_wing_ft / 100, 0.61);
		//double wd = pow((1 + rambda) / (2 * tc), 0.36);
		//double wf = pow(1 + Ve_knot / 500, 0.5);

		double S_w = A_wing * 2;				// 윗면적이 아닌, 노출 면적.
		double W_fw = W_bat;					// 연료 무게 -> 배터리 무게로 대체.
		double A = AR;							// 가로세로비.
		double Rambda = 3.5 * pitorad;			// Sweep Back Angle.
		double rambda = 0.8;					// taper ratio. 0.8로 가정.
		double q = 0.5 * rho0 * pow(V_max, 2);	// 동압.
		double N = 3.5;							// 안전계수.
		double tc = 0.12044;					// thick chord.
		double w_dg = W_gross;

		double W_wing = 0.036 *
			pow(S_w, 0.758) *
			pow(W_fw, 0.0035) *
			pow(A / cos(Rambda), 0.6) *
			pow(q, 0.006) *
			pow(rambda, 0.04) *
			pow(100 * tc / cos(Rambda), -0.3) *
			pow(N * w_dg, 0.49);
		printf("W_wing = %lf [kg]\n", W_wing);

		double L = 6;		// 동체 길이.
		double L_ft = L * mtoft;
		double W = 1.5;		// 동체 최대 폭.
		double W_ft = W * mtoft;
		double D = 1.5;		// 동체 최대 높이.
		double D_ft = D * mtoft;

		double fa = pow(W_gross_lb * N / pow(10, 5), 0.286);
		double fb = pow(L_ft / 10, 0.857);
		double fc = (W_ft + D_ft) / 10;
		double fd = pow(Ve_knot / 100, 0.338);
		double W_fuse = 200 * pow(fa * fb * fc * fd, 1.1) * lbtokg;
		printf("W_fuse = %lf [kg]\n", W_fuse);

		double tv = 0.1; // tail thickness. NACA0010
		double A_vt_ft = A_vt * pow(mtoft, 2);
		double vt_span_ft = vt_span * mtoft;

		double va = pow(W_gross_lb * N / pow(10, 5), 0.87);
		double vb = pow(A_vt_ft * 2 / 100, 1.2);
		double vc = pow(vt_span_ft / tv, 0.5);

		double W_vt = 98.5 * pow(va * vb * 0.289 * vc, 0.458) * lbtokg;
		printf("W_vt = % lf[kg] \n", W_vt);

		double W_seat = 15 * 5;
		printf("W_seat = %lf [kg] \n", W_seat);

		double W_lg = W_gross * 0.03;
		printf("W_lg = %lf [kg] \n", W_lg);

		double W_rod = 137.497 * 2;
		printf("W_rod = %lf [kg]\n", W_rod);

		double W_prop = 30 * 7;
		printf("W_prop = %lf [kg] \n", W_prop);

		double W_etc = W_gross * 0.08;
		printf("W_etc = %lf [kg] \n", W_etc);

		double W_empty = W_bat + W_wing + W_fuse + W_vt + W_seat + W_lg + W_motor + W_rod + W_prop + W_etc;
		printf("W_empty = %lf[kg] \n\n", W_empty);

		payload_calc = W_gross - W_empty;

		printf("W_gross = %lf [kg]\n", W_gross);
		printf("payload_calc = %lf [kg]\n\n\n\n\n", payload_calc);


		if (payload_calc > W_payload_req)
			W_gross = W_gross * 0.5;
		else
			W_gross = W_gross * 1.5;
	}

	return 0;
}