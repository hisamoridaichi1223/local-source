/* MPS法コード "mps.c" */

#define _CRT_SECURE_NO_WARNINGS //セキュリティが不十分であるとのエラーを表示しないため

#include <stdio.h>              //特定の関数を使えるようにするため
#include <sys/types.h>          //特定の関数を使えるようにするため
#include <math.h>               //数学的関数を使えるようにするため
#include <string.h>             //特定の関数を使えるようにするため
#include <time.h>               //時間的関数を使えるようにするため
#include <stdlib.h>             //特定の関数を使えるようにするため

// 関数全体で定義された値
#define YES 0
#define NO 1

#define NN 50000
#define DIM2 3
#define DIM 3
#define FLUID 5
//読み込みファイル
#define IN_PARAM "mps.data"  //物性値の読み込み
#define IN_PROF "mps.grid"   //粒子の初期配置読み込み
//出力ファイル
#define OUT_PROF "try1-mps.MGF" //粒子位置の出力
#define OUT_PROF4 "pressure.MGF"  //圧力の出力
#define OUT_PROF9 "wave hight.MGF" //波の高さの出力
#define OUT_PROF10 "velocity.MGF"  //速度の出力
#define OUT_PROF11 "The amount of water over pass break water.MGF" //防波堤を越えた液体の量の出力

#define NEIGHBOR 500
#define GHOST -1 //ゴースト粒子は-1
#define DT_TOO_SMALL 0.0001 //時間刻み幅？
#define INCREASE_DT 1.1　　
#define PI 3.141592  //円周率の値の定義
#define FG 10
#define flength 13.0
#define PRESSmax 2000.0  //圧力の最大値
#define PRESSmin 0.0     //圧力の最小値
#define TEMPRERATUREmax 100.0  //温度の最大値
#define TEMPRERATUREmin 0.0    //温度の最小値
#define DTMAX 0.01
#define WG 10
#define DIA 0.02  //粒径、粒子間距離
//関数の定義
void cal_n();
double weight();  //重さの関数
double cal_dt();
void testing();
void cal_buoyancy(); //浮力の関数
void cal_hyomen();
void cal_convection(); //対流の関数
void cal_confirm();
void cal_mem();
void cal_convection2();
void set_source();
void set_matrix();
void set_bcon();
void check_bcon();
void incom_decomp();
int pcg_solver();
void solver_ll();
void rev_pressgrad();
void rev_pressgrad2();
void rev_pressgrad3();
void rev_pressgrad4();
void rev_convection();
void reset_neigh();
void set_neigh();
void file1input();
void file2input();
void beoutput();
void fileoutput();
void fileoutput2();
void timeovercheck();
FILE *filecheckopen();
void set_vector_zero();
void copy_vectors();
void add_vectors();
void sub_vectors();
void mul_matrix_vector();
void mul_vectors();
void collision();
int checkconv();
int outputcheck();
int outputcount = 0;
int outputcount2 = 0;
int outputcount3 = 0;
int outputcount4 = 0;
void exit_job();
void set_pmin();
int dt_too_small();
double rho2pg();
double rho2sm();
double rho2harmonic();
void jushin();
void jushin2();
double kanseiMoment();
void jushinIdou();
double rotation();
void displacement();
void volatilization();
void testing2();
void cal_temp();
void fileoutput3();
void fileoutput4();
void fileoutput5();
void set_neigh_temp();
void fileoutput6();
void fileoutput8();
void fileoutput9();
void fileoutput10();
void fileoutput11();
void fileoutput12();
void fileoutput13();
void fileoutput14();
void fileoutput15();
double set_abc();
double set_ac();
void set_abc2();	
void Inverse();
void Inverse2d();	
void set_Vir();
void set_Vir2d();	
void set_neigh_wall();
void set_neigh_wall2();
void set_neigh_wall2d();
void penetrationcut();
void wall_edge();
void wall_move();
void supress_explosion();

//メインループ関数
void main()
{
	int fluid; //液体の種類の番号
	double rho[FLUID]; //密度
	double alpha[FLUID]; /* compressibility */
	double grav; //重力の値
	double kvis[FLUID];
	double c1[FLUID], c2[FLUID];
	static double p_surface[NN], p_wettability[NN];
	static double prwall_total1[DIM][NN], pr_total1[DIM][NN];
	static double prwall_total2[DIM][NN], pr_total2[DIM][NN];
	double tmax, dtmx, dtra, epsp; /* computation parameters */
	double disa, rerp, rericcg; /* kernel parameters */
	int itmx, imxp1, imnp1; /* computation parameters */
	char outputtype[3]; /* output time; I : iteration, T : time */
	int interval; //出力の頻度
	double output_dt, output_time; /* output interval, output time */
	double time_sim; /* simulation time [sec] */
	static double p[NN];//p_surface[NN],p_wettability[NN]
	static double xn[DIM][NN], rd[DIM][NN], vn[DIM][NN], dv[DIM][NN], dp[NN], source[NN], lv[DIM][NN];
	static double n[NN], n_fluid[FLUID][NN], n2[NN];
	static double prwall_total[DIM][NN], pr_total[DIM][NN];
	static int typep[NN]; /* particle type */
	int nump, numf, numw, num_fluid, numg; /* 総粒子数 */
	double rep, rep2, reiccg, reiccg2; /* re of the kernel */
	double n0p, n0iccg; /* fixed particle number density */
	int n0p_particle; /* particle to evaluate n0p in the initial configuration */
	double dt,dt2;
	double lambda; /* int (w r^2 dv) / int (w dv) */
	int it;
	static int neigh[NN][NEIGHBOR]; /* neighboring particle numbers storage */
	static int neigh_iccg[NN][NEIGHBOR];
	static double poiss[NN][NEIGHBOR], ic[NN][NEIGHBOR];
	static int bcon[NN]; /* -1 : ignore (wall,symmetry), 0 : to be solved, 1 : value fixed */
	static int check[NN];
	static int start_i, next_i;
	static int count1;
	double r2lim;
	static int group[NN];
	static int fgroup[NN];
	int i, j, k, l, kk, a, lll;
	FILE *fp1, *fp2, *fpo, *fpo4, *fpo9; //ファイルの設定
	double minlen, dirichlet, col_rat, dirichletcut;
	int num_dirichlet = 0;
	double ra;
	int wall_type, wall_p_type;
	static double pmin[NN];
	double lam_rat;
	//	static double rg[DIM][FG],rg2[DIM][FG],rgn[DIM][FG];
	double rot[FG], ang[FG];
	double KM;
	double tempn[NN], ltemp[NN];
	double cond[FLUID], cp[FLUID];
	double xn_vir[DIM][NN];							//追加（：仮想粒子）
	double w1[DIM][WG], w2[DIM][WG], w3[DIM][WG];
	double w4[DIM][WG], w5[DIM][WG], w6[DIM][WG];	//追加（全て壁の数を上限100と仮定）
	double Nor_vec[DIM][WG];						//追加(：法線ベクトル n)
	double Inv[3][3][WG];								//追加(：逆行列の各値)
	double xn_wall[DIM][NN][NEIGHBOR];	//追加
	int neigh_wall[NN];
	int bou_flag[NN][NEIGHBOR];
	double k2[WG], k3[WG];
	double vmax;
	double wall_pene[WG][2][3];
	double vmaxcut, vt, pmax;
	int overpass[NN], waterover[2], oilover[2];//超えたら1で超えてなかったら0
	double acceleration, maxv, start, finish, deceleration;

	dt = 0;
	a = 0;

/*      乱数      */
srand((unsigned)time(NULL));

/*ファイルの設定*/

fprint(stderr, "%%%% mps %%%%\n");

/*入力ファイルを設定し開く*/

	fp1 = filecheckopen(1, IN_PARAM, "r");
	fp2 = filecheckopen(2, IN_PROF, "r");
/*出力ファイルを設定し開く*/
	fpo = filecheckopen(3, OUT_PROF, "w");
	fpo4 = filecheckopen(6, OUT_PROF4, "w");
	fpo9 = filecheckopen(11, OUT_PROF9, "w");
        fpo11 = filecheckopen(13, OUT_PROF11, "w");

/* 入力ファイルからのデータの読み込み */

	file1input(fp1, &fluid, &wall_type, &wall_p_type, rho, alpha, &grav,
		&tmax, &itmx, &dtmx, &dtra, &epsp, &imxp1, &imnp1, outputtype, &interval, &output_dt, &disa, &rerp, &rericcg,
		&minlen, &dirichlet, &dirichletcut, &col_rat, &lam_rat,
		&n0p, &n0iccg, &n0p_particle, kvis, c1, c2, cp, cond, &vmaxcut, &pmax, &acceleration, &maxv, &start, &finish, &deceleration);

	file2input(fp2, &time_sim, &nump, &numf, &numw, &typep, xn, vn, p, n, &fgroup, tempn, w1, w2, w3, Nor_vec, wall_pene, &numg);


	/* 計算パラメータの設定 */

	rep = disa*rerp; /* real distance [m] of re(PND) */
	rep2 = rep*rep;
	reiccg = disa*rericcg; /* real distance [m] of re(ICCG) */
	reiccg2 = reiccg*reiccg;

        /* 液体粒子の数をカウント */
	vmax = 0;
	num_fluid = 0;
	vt = 0;

	for (i = 0; i < nump; i++) {
		if (typep[i] == 0 || typep[i] == 4 || typep[i] == GHOST)num_fluid++;
	}
        /* 防波堤を越えた粒子数に0を設定 */
	for (i = 0; i < nump; i++) {
		overpass[i] = 0;
	}
	waterover[0] = 0;
	oilover[0] = 0;
	waterover[1] = 0;
	oilover[1] = 0;

	/* If n0p=0.0 in the input data, n0p is given as n[n0p_particle] */
	/* "n0p_particle" is given in the input data */

	if (n0p == 0.0) {
		reset_neigh(nump, neigh);
		reset_neigh(nump, neigh_iccg);
		for (i = 0; i<nump; i++) {
			if (typep[i] == GHOST) continue;
			set_neigh(i, xn, nump, typep, neigh, rep2, neigh_iccg, reiccg2, wall_type);
			set_neigh_wall2d(i, Nor_vec, Inv, w1, xn_vir, xn_wall, neigh_wall, numw, xn, rep, bou_flag, a, dt, k2, k3, disa);
		}
		cal_n(i, xn, n, n2, n_fluid, neigh, typep, rep, fluid, numw, xn_wall, neigh_wall);
		n0p = n[n0p_particle];
		cal_n(i, xn, n, n2, n_fluid, neigh, typep, rep, fluid, numw, xn_wall, neigh_wall);
		n0iccg = n[n0p_particle];
		fprintf(stderr, "calculated n0p, n0iccg = %e, %e\n", n0p, n0iccg);
	}

	/* Particles collide when r < r2lim */

	r2lim = disa*disa*disa*minlen*minlen*minlen;//三次元なら三つ

	/* pi * ra^2 = disa^2 */

	ra = disa / 1.7724539;

	lambda = lam_rat*(reiccg2*reiccg2 / 12.0 - reiccg*ra*ra*ra / 3.0 + ra*ra*ra*ra / 4.0)
		/ (reiccg2 / 2.0 - reiccg*ra + ra*ra / 2.0);

	fprintf(stderr, "rep,reiccg=%e,%e\n", rep, reiccg);
	fprintf(stderr, "r2lim=%e\n", r2lim);
	fprintf(stderr, "lambda=%e\n", lambda);

	/* Initla dt setting */
	

	dt = dtmx;
	dt2=dt;
	a = 1;


	/* Initial setting of neighboring tables */

	reset_neigh(nump, neigh);
	reset_neigh(nump, neigh_iccg);
	for (i = 0; i<nump; i++) {
		if (typep[i] == GHOST) continue;
		set_neigh(i, xn, nump, typep, neigh, rep2, neigh_iccg, reiccg2, wall_type);
		//set_neigh_wall2d(i, Nor_vec, Inv, w1, xn_vir, xn_wall, neigh_wall, numw, xn, rep, bou_flag, a, dt, k2, k3, disa);
		cal_n(i, xn, n, n2, n_fluid, neigh, typep, rep, fluid, numw, xn_wall, neigh_wall);
	}

	/* Initial setting of variables */





	/*  particle number density */

	//	for (i = 0; i<nump; i++) {
	//		if (typep[i] == GHOST || typep[i] == wall_type) continue;
	//cal_n(i,xn,n,n_fluid,neigh,typep,rep,fluid);
	//	}

	/*  boundary codition setting */

	for (i = 0; i<nump; i++) {
		set_bcon(i, typep, bcon, n, n0p, dirichlet, wall_type);
	}


	/*  source term */

	//	num_dirichlet = 0;
	for (i = 0; i<nump; i++) {
		if (typep[i] == GHOST || typep[i] == wall_type) continue;
		if (bcon[i] == 0)set_source(i, typep, source, n, n0p, rho, dt);
		else if (bcon[i] == 1) {
			//source[i]=0.0;
			set_source(i, typep, source, n, n0p, rho, dt);
		}
	}
	for (i = 0; i<nump; i++) {
		if (p[i]>pmax) {
			p[i] = pmax;
		}
	}

	check_bcon(nump, typep, bcon, neigh, check);



	/* Next output timing */

	if (outputtype[0] == 'T')output_time = time_sim + output_dt;



	/* Initial particles' profile is output to %%%% file */
	beoutput(fpo);
	beoutput(fpo4);
	a = 1;
	fileoutput(fpo, time_sim, nump, typep, xn, vn, p, n, bcon, a, disa,fgroup);
	fileoutput4(fpo4, time_sim, nump, typep, xn, vn, p, n, bcon, a, disa, num_fluid);
	fileoutput9(fpo9, nump, typep, time_sim, xn);
	fileoutput11(fpo11, nump, typep, time_sim, xn, overpass, waterover, oilover);


/*メインループ*/

	for (it = 1; it <= itmx; it++) {

		/*  刻み幅、時間経過 */
		
		dt2 = cal_dt(typep, vn, xn, nump, dtra, disa, it, dtmx, wall_type, wall_p_type, dt2, INCREASE_DT, &vmax, &vmaxcut, n, &dirichletcut);
		if(time_sim<finish)dt=dt2/2;
		else dt=dt2;


		time_sim = time_sim + dt;

		
		/* 刻み幅が小さすぎると解析が止まる設定 */

		if (dt_too_small(dt2, dtmx, DT_TOO_SMALL) == YES) {
			fprintf(stderr, "**** The time step falls to a too small value, calculaion ends.  dt = %f\n", dt);
			exit_job();  
		}


		/*　移動壁の関数　*/
		wall_move(nump,fgroup,xn,vn,dt,it,&vt,time_sim, acceleration, maxv, start, finish, deceleration);



		/*  setting of neighboring tables */

		reset_neigh(nump, neigh);
		reset_neigh(nump, neigh_iccg);
		for (i = 0; i<nump; i++) {
			if (typep[i] == GHOST) continue;
			set_neigh(i, xn, nump, typep, neigh, rep2, neigh_iccg, reiccg2, wall_type);
			



		/** seach surface particle   ********/
		/**         and              ********/
		/** wall particle has been removed **/

		for (i = 0; i<nump; i++) {
			if (typep[i] == GHOST || typep[i] == wall_type || typep[i] == wall_p_type) continue;
			cal_n(i, xn, n, n2, n_fluid, neigh, typep, rep, fluid, numw, xn_wall, neigh_wall);
			set_bcon(i, typep, bcon, n, n0p, dirichlet, wall_type);
		}

		if (it == 179) {
			fprintf(stderr, "#116 ");
		}


		//速度項のラプラシアンの計算(粘性項)
		set_vector_zero(nump, lv[0]);
		set_vector_zero(nump, lv[1]);
		set_vector_zero(nump, lv[2]);
		for (i = 0; i<nump; i++) {
			//		if ( typep[i]==GHOST || typep[i]==wall_type || typep[i]==wall_p_type ) continue;
			if (typep[i] == GHOST) continue;
			testing(i, typep, neigh, n, rho, xn, vn, lv, rep2, rep, n0p, wall_type, lambda, xn_wall, neigh_wall);
		}

		//表面張力・濡れ性ポテンシャル関数の決定
		for (i = 0; i<nump; i++) {
		if (typep[i] == GHOST || typep[i] == wall_type || typep[i] == wall_p_type) continue;
		//cal_hyomen(i, disa, rep, wall_type, wall_p_type, neigh, typep, dt, xn, p_surface, p_wettability, prwall_total, pr_total, bcon,c1,c2);
		}
		






		/*******************/
		/*	仮の速度を決定 */
		/*******************/

		/*  buoyancy force */
		for (i = 0; i<nump; i++) {
			if (typep[i] == GHOST) continue;
			
			cal_buoyancy(i, lv, vn, grav, dt, prwall_total, pr_total, typep, rho, fgroup, kvis, c1, c2,disa);
		}


		/*  convection (particle motion) */

		for (i = 0; i<nump; i++) {
			if (typep[i] == GHOST || typep[i] == wall_type || typep[i] == wall_p_type) continue;
			cal_convection(i, xn, vn, dt);
		}

		for (i = 0; i<nump; i++) {
			penetrationcut(i, xn, typep, wall_pene, numg);//ゴースト空間内の粒子を消す
		}

		collision(nump, typep, bcon, neigh, xn, vn, rho, r2lim, dt, rep2, col_rat, wall_type, wall_p_type);




		/*  setting of neighboring tables */

		reset_neigh(nump, neigh);
		reset_neigh(nump, neigh_iccg);
		for (i = 0; i<nump; i++) {
			if (typep[i] == GHOST) continue;
			set_neigh(i, xn, nump, typep, neigh, rep2, neigh_iccg, reiccg2, wall_type);
			//				set_neigh_wall2d(i, Nor_vec, Inv, w1, xn_vir, xn_wall, neigh_wall, numw, xn, rep, bou_flag, a, dt, k2, k3, disa);
		}


		/*                       */
		/*  pressure correction  */
		/*        ICCG           */

		/*  particle number density */

		for (i = 0; i<nump; i++) {
			if (typep[i] == GHOST || typep[i] == wall_type) continue;
			cal_n(i, xn, n, n2, n_fluid, neigh, typep, rep, fluid, numw, xn_wall, neigh_wall);
		}

		/*  boundary codition setting */

		for (i = 0; i<nump; i++) {
			set_bcon(i, typep, bcon, n, n0p, dirichlet, wall_type);
		}

		/*  source term */

		num_dirichlet = 0;
		for (i = 0; i<nump; i++) {
			if (typep[i] == GHOST || typep[i] == wall_type) continue;
			if (bcon[i] == 0)set_source(i, typep, source, n, n0p, rho, dt);
			else if (bcon[i] == 1) {
				set_source(i, typep, source, n, n0p, rho, dt);
			}
		}

		//		check_bcon(nump, typep, bcon, neigh, check);

		/*  coefficient matrix */

		for (i = 0; i<nump; i++) {
			set_matrix(i, typep, bcon, neigh_iccg, poiss, xn, n0iccg, reiccg, reiccg2, lambda, rho, alpha, dt);
		}


		/*  incomplete Choleskey decomposition */

		incom_decomp(nump, bcon, neigh_iccg, poiss, ic);

		/*  preconditioned conjugate gradient solver */

		set_vector_zero(nump, dp);

		/* ICCG iteration */

		pcg_solver(nump, bcon, neigh_iccg, poiss, dp, source, ic, epsp, imnp1, imxp1, dt);

		/*  minus pressure -> zero pressure */

		for (i = 0; i<nump; i++) {
			if (dp[i]<0.0)dp[i] = 0.0;
		}


		/*  maximum pressure table generation */

		for (i = 0; i<nump; i++) {
			if (bcon[i] == -1 || typep[i] == wall_p_type) continue;
			set_pmin(i, typep, neigh_iccg, dp, pmin, wall_type);
		}
		//	supress_explosion(bcon,dp,neigh_iccg,nump);

		/*  velocity correction, pressure gradient term */

		set_vector_zero(nump, dv[0]);
		set_vector_zero(nump, dv[1]);
		set_vector_zero(nump, dv[2]);


		for (i = 0; i<nump; i++) {
			if (bcon[i] == -1 || typep[i] == wall_p_type) continue;
			
			rev_pressgrad2(i, typep, neigh_iccg, n, rho, xn, dv, dp, dt, rep2, rep, n0p, wall_type, pmin, numw, xn_wall, neigh_wall, disa);
		}


		/*
		add_vectors(nump,vn[0],vn[0],dv[0],1);
		add_vectors(nump,vn[1],vn[1],dv[1],1);
		add_vectors(nump,vn[2],vn[2],dv[2],1);
		*/
		copy_vectors(nump, dp, p);




		/*  correction of particle motion */

		for (i = 0; i<nump; i++) {
			if (bcon[i] == -1 || typep[i] == wall_p_type) continue;
			cal_convection(i, xn, dv, dt);
		}

		for (i = 0; i<nump; i++) {
			vn[0][i] = vn[0][i] + dv[0][i];
			vn[1][i] = vn[1][i] + dv[1][i];
				vn[2][i] = vn[2][i] + dv[2][i];//三次元ならいる
		}



		/*  出力設定  */

		if (outputcheck(outputtype, time_sim, &output_time, output_dt, it, interval, a) == YES)
		{
			a = a + 1;
			fileoutput(fpo, time_sim, nump, typep, xn, vn, p, n, bcon, a, disa,fgroup);
			fileoutput4(fpo4, time_sim, nump, typep, xn, vn, p, n, bcon, a, disa, num_fluid);
			fileoutput9(fpo9, nump, typep, time_sim, xn);
			fileoutput11(fpo11, nump, typep, time_sim, xn, overpass, waterover, oilover);
		}



		/*  copy xn -> x */


		timeovercheck(time_sim, tmax, a);


	} /* end of time loop */

	fclose(fpo2);

}