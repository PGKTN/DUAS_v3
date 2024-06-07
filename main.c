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
	
	double W_gross = 3000;
	double W_payload_req = 550;
	double payload_calc = 0.0;

	int i = 0;

	while (abs(W_payload_req - payload_calc) > 10e-6)
	{
		printf("----------------%d번째----------------\n", i);
		i++;


		double wing_span = 9.402822;
		double wing_root = 1.8;
		double wing_tip = 1.2;
		double wing_chord = (wing_root + wing_tip) / 2;
		double AR = wing_span / wing_root;
		printf("AR = %lf \n", AR);
		double A_wing = wing_span * wing_chord;

		double WL = W_gross / A_wing;

		double r_disk = 1.5;
		double A_disk = Pi * pow(r_disk, 2);
		double N_rotor = 4;
		double N_crotor = 3;









		double T_hover = W_gross * g * 1.2;

		double P_hover = sqrt(pow(T_hover, 3) / (2 * rho0 * A_disk * 4));
		double pHOVER = P_hover / W_gross;
		printf("P_hover = %lf [W]\n", P_hover);
		double V_max = 55.55555556;

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







		////////////////////////
		/// Mission Analysis ///
		////////////////////////




		double range_cruise = 50 * 1000; // [m]
		double t_cruise = range_cruise / V_max;
		double rho_cruise = Get_rho(rho0, 600);
		double W_gross_lb = W_gross * kgtolb;

		double fe_fuse = 1.7 * pow(((W_gross_lb) / 1000), 2 / 3);		// lb -> kg
		double D_fuselage = fe_fuse * 0.5 * rho_cruise * pow(V_max, 2);

		double fe_hub = 0.85 * pow(((W_gross_lb) / 1000), 2 / 3);		// lb -> kg
		double D_hub = fe_hub * 0.5 * rho_cruise * pow(V_max, 2);

		double alpha = 4.3 * Pi / 180; // degree to radian
		double Cl = 0.9599;
		double L_wing = 0.5 * rho_cruise * pow(V_max, 2) * A_wing * Cl;
		//printf("L_wing = %lf [N]\n", L_wing);
		double Cd = 0.00542;
		double D_wing = 0.5 * rho_cruise * pow(V_max, 2) * A_wing * Cd;
		//printf("D_wing = %lf [N]\n", D_wing);

		double F_x = D_hub * 4 + D_fuselage + D_wing * cos(alpha) + L_wing * sin(alpha);
		double T1 = W_gross * g - L_wing * cos(alpha) + D_wing * sin(alpha);

		while (abs(T1) > 10e-6)
		{
			if ((L_wing * cos(alpha) - D_wing * sin(alpha)) > W_gross * g)
				wing_span = wing_span * 0.5;
			else
				wing_span = wing_span * 1.5;

			A_wing = wing_span * wing_chord;
			L_wing = 0.5 * rho_cruise * pow(V_max, 2) * A_wing * Cl;
			D_wing = 0.5 * rho_cruise * pow(V_max, 2) * A_wing * Cd;
			F_x = D_hub * 4 + D_fuselage + D_wing * cos(alpha) + L_wing * sin(alpha);
			T1 = W_gross * g - L_wing * cos(alpha) + D_wing * sin(alpha);
		}

		printf("wing_span = %lf [m]\n", wing_span);
		printf("T1 = %lf [N] \n", T1);

		double P_cruise_total = F_x * V_max;
		printf("P_cruise_total = %lf [W] \n", P_cruise_total);
		double P_cruise_motor = P_cruise_total / N_crotor;
		printf("P_cruise_motor = %lf [W]\n", P_cruise_motor);
		double E_cruise = P_cruise_total * t_cruise / 3600;
		printf("E_cruise = %lf [Wh]\n\n", E_cruise);




		double T = W_gross * g;
		double T_rotor = T / N_rotor;

		double sigma = 0.2262;

		double rpm_climb = 1058;
		double omega_climb = Get_omega(rpm_climb);
		double v_tip_climb = omega_climb * r_disk;
		double rho_climb = Get_rho(rho0, 0);
		double t_climb = 3 * 60;

		double Cd0 = 0.03390;

		double P_climb_total = 0.0;
		double E_climb = 0.0;


		double vi_climb = Get_vi(T_rotor, rho_climb, A_disk);
		//printf("vi_climb = %lf [m/s]\n", vi_climb);
		double v_climb = 3.333333;
		double P_i_climb = Get_P_i(T_rotor, vi_climb + v_climb);
		//printf("P_i_climb = %lf [W]\n", P_i_climb);
		double P_0_climb = Get_P_0(sigma, Cd0, rho_climb, A_disk, v_tip_climb);
		//printf("P_0_climb = %lf [W]\n", P_0_climb);
		double P_climb_motor = P_i_climb + P_0_climb;
		printf("P_climb_motor = %lf [W]\n", P_climb_motor);
		P_climb_total = P_climb_motor * N_rotor;
		printf("P_climb_total = %lf [W]\n", P_climb_total);
		E_climb = P_climb_total * (t_climb / 3600);
		printf("E_climb = %lf [Wh]\n\n", E_climb);


		double t_hover = 0;
		double rpm_hover = 1058;
		double omega_hover = Get_omega(rpm_hover);
		double v_tip_hover = omega_hover * r_disk;
		double vi = Get_vi(T_rotor, rho_climb, A_disk);

		double P_i_hover = Get_P_i(T_rotor, vi);
		double P_0_hover = Get_P_0(sigma, Cd0, rho0, A_disk, v_tip_hover);

		double P_hover_motor = P_i_hover + P_0_hover;
		double P_hover_total = P_hover_motor * N_rotor;

		while (abs(P_hover_total - P_hover) > 1)
		{
			if (P_hover_total > P_hover)
				rpm_hover = rpm_hover * 0.5;
			else
				rpm_hover = rpm_hover * 1.5;

			omega_hover = Get_omega(rpm_hover);
			v_tip_hover = omega_hover * r_disk;
			P_i_hover = Get_P_i(T_rotor, vi);
			P_0_hover = Get_P_0(sigma, Cd0, rho0, A_disk, v_tip_hover);

			P_hover_motor = P_i_hover + P_0_hover;
			P_hover_total = P_hover_motor * N_rotor;
		}

		printf("rpm_hover = %lf\n", rpm_hover);	// 1266

		printf("P_hover_motor = %lf [W]\n", P_hover_motor); // 199401 [W]
		printf("P_hover_total = %lf [W]\n", P_hover_total);
		double E_hover = P_hover_total * t_hover / 3600;
		printf("E_hover = %lf [Wh]\n\n", E_hover);



		double rpm_descent = 1058;
		double omega_descent = Get_omega(rpm_climb);
		double v_tip_descent = omega_descent * r_disk;

		double rho_descent = Get_rho(rho0, 600);
		double t_descent = 9 * 60;
		double vi_descent = Get_vi(T_rotor, rho_descent, A_disk);
		double v_descent = -1.111111;
		double P_i_descent = Get_P_i(T_rotor, vi_descent + v_descent);
		double P_0_descent = Get_P_0(sigma, Cd0, rho_descent, A_disk, v_tip_descent);
		double P_descent_motor = P_i_descent + P_0_descent;
		printf("P_descent_motor = %lf [W]\n", P_descent_motor);
		double P_descent_total = P_descent_motor * N_rotor;
		printf("P_descent_total = %lf [W]\n", P_descent_total);
		double E_descent = P_descent_total * (t_descent / 3600);
		printf("E_descent = %lf [Wh]\n\n", E_descent);



		double E_total = E_hover + E_climb + E_cruise + E_descent;

		double density_bat = 380; // [Wh/kg]
		double eta_prop = 0.8;
		double W_bat = E_total * 1.2 / (density_bat * eta_prop);
		printf("W_bat = %lf [kg]\n", W_bat);


		double N = 3.4; // 하중계수.
		double Rambda = 3.5 * pitorad;
		double rambda = wing_tip / wing_root;
		double tc = 0.12044;
		rho_cruise = Get_rho(rho0, 600);
		double V_max_knot = V_max * mstoknot;
		double Ve_knot = 0.5 * rho_cruise * pow(V_max_knot, 2);
		W_gross_lb = W_gross * kgtolb;
		double A_wing_ft = A_wing * pow(mtoft, 2);

		double wa = pow(W_gross_lb * N / 10e5, 0.65);
		double wb = pow(AR / cos(Rambda), 0.57);
		double wc = pow(A_wing_ft / 100, 0.61);
		double wd = pow((1 + rambda) / (2 * tc), 0.36);
		double wf = pow(1 + Ve_knot / 500, 0.5);

		double W_wing = 77.6 * pow(wa * wb * wc * wd * wf, 0.993) * lbtokg;
		printf("W_wing = %lf [kg]\n", W_wing);

		double L = 6;		// 동체 길이.
		double L_ft = L * mtoft;
		double W = 1.5;		// 동체 최대 폭.
		double W_ft = W * mtoft;
		double D = 1.5;		// 동체 최대 높이.
		double D_ft = D * mtoft;

		double fa = pow(W_gross_lb * N / 10e5, 0.286);
		double fb = pow(L_ft / 10, 0.857);
		double fc = (W_ft + D_ft) / 10;
		double fd = pow(Ve_knot / 100, 0.338);
		double W_fuse = 200 * pow(fa * fb * fc * fd, 1.1) * lbtokg;
		printf("W_fuse = %lf [kg]\n", W_fuse);

		double W_seat = 15 * 5;
		printf("W_seat = %lf [kg] \n", W_seat);

		double W_lg = W_gross * 0.03;
		printf("W_lg = %lf [kg] \n", W_lg);

		double W_motor = 18.7 * 7;
		printf("W_motor = %lf [kg] \n", W_motor);

		double W_etc = W_gross * 0.10;
		printf("W_etc = %lf [kg] \n", W_etc);

		double W_empty = W_bat + W_wing + W_fuse + W_seat + W_lg + W_motor + W_etc;
		printf("W_empty = %lf[kg] \n\n", W_empty);

		payload_calc = W_gross - W_empty;

		printf("W_gross = %lf [kg]\n", W_gross);
		printf("payload_calc = %lf [kg]\n", payload_calc);


		if (payload_calc > W_payload_req)
			W_gross = W_gross * 0.5;
		else
			W_gross = W_gross * 1.5;
	}



	return 0;
}