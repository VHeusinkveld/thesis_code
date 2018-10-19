struct sDiag dia; 		// Diagnostics

struct sDiag {
	double Ekin;			// Total kinetic energy
	double EkinOld;			// Track changes in kin energy 
	double WdoneOld;		// Track changes in work done 
	double rotVol;			// Track real rotor volume
};

/* Sanity checks */
event sanity (t += 1.){
	scalar ekin[]; 		// Kinetic energy
	double tempVol = 0.;
	double tempEkin = 0.;	
	double max_vel = 0.;
	
	foreach(reduction(+:tempVol) reduction(+:tempEkin) reduction(max:max_vel)) {
		#if dimension == 2			
			ekin[] = 0.5*rho[]*cube(Delta)*(sq(u.x[]) + sq(u.y[])));
			tempVol += sq(Delta)*fan[];
		#elif dimension == 3
			ekin[] = 0.5*rho[]*cube(Delta)*(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
			tempVol += cube(Delta)*fan[];
		#endif
		tempEkin += ekin[];
		max_vel = max(max_vel, sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[])));
	}	

	dia.rotVol = 1.*tempVol;
	dia.Ekin = 1.*tempEkin;
	
	if(t <= 10.) {
		printf("eps=%g, rotv=%g, maxv=%g\n", eps, rot.cu, max_vel);
	}

	if(fabs(dia.rotVol/rot.V - 1) > 0.1){
		printf("V=%g, Vr=%g, ",rot.V, dia.rotVol);
	}	
	printf("Energy: Ek=%g, W=%g, Ek/W=%g, dEk/dW=%g\n", 
		dia.Ekin, rot.Work, dia.Ekin/rot.Work, 
		(dia.Ekin-dia.EkinOld)/(rot.Work-dia.WdoneOld));
	
	dia.EkinOld = 1.*dia.Ekin;
	dia.WdoneOld = 1.*rot.Work;
}
