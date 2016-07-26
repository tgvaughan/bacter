package bacter.devutils;

import beast.core.BEASTObject;

public class Buffer extends BEASTObject {
	
	public double [][][] cfTransitionProbs;
	
	public void iniDimensions(int x, int y, int z){
		cfTransitionProbs = new double[x][y][z];
		System.err.println("Initialize dimemnsions "+ x + " "+y+ " "+ z);
	}

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}

}
