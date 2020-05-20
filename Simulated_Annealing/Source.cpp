#include<iostream>
#include<random>

#include"Simulated_Annealing.h"
using namespace std;


int main()
{
	initialize_rand();
	Sim_Annealing FloorPlanAnneal;
	int Plot_N = 1;
	FloorPlanAnneal.InputfileRead();
	FloorPlanAnneal.Init_Polish_Exp();
	float  cost;
	FloorPlanAnneal.T_final = 0.5f;
	FloorPlanAnneal.last_accepted_cost = FloorPlanAnneal.Calculate_Anneal_Cost();
	FloorPlanAnneal.Polish_Expressionbest = FloorPlanAnneal.Polish_Exp;
	FloorPlanAnneal.Anneal_move2();
	cost = FloorPlanAnneal.Calculate_Anneal_Cost();
	FloorPlanAnneal.T = (-abs(FloorPlanAnneal.last_accepted_cost - cost)) / log(0.97);
	printf("Initial Temperature: %f\n", FloorPlanAnneal.T);
	printf("\n");
	printf("Initial Polish Expression: \n");
	FloorPlanAnneal.printPolish_Exp();
	printf("\n");
	if (cost < FloorPlanAnneal.last_accepted_cost)
	{
		FloorPlanAnneal.Polish_Expressionbest = FloorPlanAnneal.Polish_Exp;
		FloorPlanAnneal.last_accepted_cost = cost;
	}
	else
	{
		FloorPlanAnneal.Polish_Exp = FloorPlanAnneal.Polish_Expressionbest;
		FloorPlanAnneal.Calculate_Anneal_Cost();
	}
	FloorPlanAnneal.numMoves = FloorPlanAnneal.Polish_Exp.N * 2;
	printf("Initial cost:%f\n", FloorPlanAnneal.Polish_Exp.final_node->h*FloorPlanAnneal.Polish_Exp.final_node->w);
	printf("\n");
	FloorPlanAnneal.T_dec_ratio = 0.85f;

	FloorPlanAnneal.FloorPLot(&Plot_N);
	Plot_N++;
	FloorPlanAnneal.Anneal_Schedule();
	FloorPlanAnneal.FloorPLot(&Plot_N);

	printf("Final Polish Expression: \n");
	FloorPlanAnneal.printPolish_Exp();
	FloorPlanAnneal.write_table();
	printf("\n");
	printf("Final cost : %f\n", FloorPlanAnneal.last_accepted_cost = FloorPlanAnneal.Polish_Exp.final_node->h*FloorPlanAnneal.Polish_Exp.final_node->w);
	printf("\n");
	printf("COMPLETE!\n");
	return 1;
}