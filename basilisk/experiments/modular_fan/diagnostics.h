struct sDiag dia; 		// Diagnostics

struct sDiag {
	double Ekin;			// Total kinetic energy
	double EkinOld;			// Track changes in kin energy 
	double WdoneOld;		// Track changes in work done 
	double rotVol;			// Track real rotor volume
};

/* Sanity checks */
event sanity (t += 1){
	
	scalar ekin[]; 		// Kinetic energy
	double tempVol = 0.;
	double tempEkin = 0.;	
	
	foreach(reduction(+:tempVol) reduction(+:tempEkin)) {
		ekin[] = 0.5*rho[]*sq(Delta)*(sq(u.x[]) + sq(u.y[]));
		tempEkin += ekin[];
		#if dimension > 1			
			tempVol += sq(Delta)*fan[];
		#elif dimension > 2
			tempVol = cube(Delta)*fan[];
			bf += u.y[]*b[]*sq(Delta);
		#endif
	}

	dia.rotVol = 1.*tempVol;
	dia.Ekin = 1.*tempEkin;
	

	if(fabs(dia.rotVol/rot.V - 1) > 0.1){
		printf("V=%g, Vr=%g, ",rot.V, dia.rotVol);
	}	
	printf("Energy: Ek=%g, W=%g, Ek/W=%g, dEk/dW=%g\n", 
		dia.Ekin, dia.Wdone, dia.Ekin/rot.Work, 
		(dia.Ekin-dia.EkinOld)/(rot.Work-dia.WdoneOld));
	
	dia.EkinOld = 1.*dia.Ekin;
	dia.WdoneOld = 1.*rot.Work;
}