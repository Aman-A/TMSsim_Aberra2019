: $Id: xtra.mod,v 1.4 2014/08/18 23:15:25 ted Exp ted $


NEURON {
	SUFFIX xtra
	RANGE es : (es = max amplitude of the potential)
	RANGE Ex,Ey,Ez
	RANGE Dx,Dy,Dz
	RANGE x, y, z, type, order
	GLOBAL stim : (stim = normalized waveform)
	POINTER ex 
}

PARAMETER {
	: default transfer resistance between stim electrodes and axon
	:rx = 1 (megohm) : mV/nA
	es = 0 (mV)
	x = 0 (1) : spatial coords
	y = 0 (1)
	z = 0 (1)
	Ex = 0 (1) : E-field components
	Ey = 0 (1)
	Ez = 0 (1)
	Dx = 0 (1) : unit vector of orientation of compartment (central-diff)
	Dy = 0 (1)
	Dz = 0 (1)
	type = 0 (1) : numbering system for morphological category of section - unassigned is 0
	order = 0 (1) : order of branch/collateral. 
}

ASSIGNED {
	v (millivolts)
	ex (millivolts)
	stim (unitless) 	
	:im (milliamp/cm2)
	:er (microvolts)
	area (micron2)
}

INITIAL {
	ex = stim*es
	:er = (10)*rx*im*area
}


BEFORE BREAKPOINT { : before each cy' = f(y,t) setup
  ex = stim*es
}

: AFTER SOLVE { : after each solution step
:   er = (10)*rx*im*area
: }

