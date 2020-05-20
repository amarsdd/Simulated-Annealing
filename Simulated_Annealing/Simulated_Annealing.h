#ifndef _SIMULATED_ANNEALING_H_
#define _SIMULATED_ANNEALING_H_
#include<iostream>
#include<random>
#include<chrono>
#include<vector>
#include<list>
#include<algorithm>
#include<string>
#include"EasyBMP.h"
#include"EasyBMP_Geometry.h"

using namespace std;
typedef std::chrono::high_resolution_clock myclock;

#define interpol_num	5

void initialize_rand();

class BdnCurve{
public:
	float w;
	float h;
	int BdnCurveID;
};

class node;
class Module{
public:
	int Original_Module_ID;
	int inPolish_ID;
	float A_;
	float r_;
	float s_;
	int S;
	int FinalBdnCurveIdx_tmp;
	int FinalBdnCurveIdx;
	float location[2];
	float center[2];
	node* parentnode;
	vector<BdnCurve> BdnCurveIdx;
};


class Operator{
public:
	int Type;
	int inPolish_ID;
	list<node> combineBdnCurve;
};

class node {
public:
	void initialize();
	float w;
	float h;
	node *next;
	list< node>::iterator parentnode;
	int operator_type;
	int node_type;
	int type_flag;
	int head_flag;
	Module *Module1;
	Module *Module2;
	int FinalBdnCurveIdx1;
	int FinalBdnCurveIdx2;
	Operator *Operator1;
	Operator *Operator2;
	int final_node1;
	int final_node2;
	float location[2];
};
class Operator_chain{
public:
	int inPolish_ID;
	int Operator_Serial;
	int Module_Serial;
	int inPolish_ID_R;
	int Operator_Serial_R;
	int Module_Serial_R;
};

class Polish_Expression{
public:
	vector<Module> ModuleSeq;
	vector<Operator> OperatorSeq;
	list<Operator_chain> Operator_chain;
	vector<int> Full_Polish_Exp;
	list< node>::iterator final_node;
	int N;
	float p;
	float q;
};

class FloorplanData{
public:
	float t;
	float area;
	float p;
};

class Sim_Annealing{
public:
	Sim_Annealing(){};
	float T_final;
	float T;
	int numMoves;
	float last_accepted_cost;
	float T_dec_ratio;
	float aspect_ratio_final;
	vector<FloorplanData> table;
	int num_move1 = 0, num_move2 = 0, num_move3 = 0;
	int tmp_;
	int temp_;
	Polish_Expression Polish_Exp;
	Polish_Expression Polish_Expressionbest;
	bool InputfileRead();
	float get_data(string line, string::size_type pos);
	void Init_Polish_Exp();
	void combineBdnCurvesType1(node* temp, int seq_id1, int seq_id2, int i, int operatorflag);
	void combineBdnCurvesType2(node* temp, int seq_id1, int seq_id2, int i, int operatorflag);
	void combineBdnCurvesType3(node* temp, int seq_id1, int seq_id2, int i, int operatorflag);
	void combineBdnCurvesType4(node* temp, int seq_id1, int seq_id2, int i, int operatorflag);
	float getAreaCost();
	void BdnCurveIdxCal(Module *Module);
	void Remove_Redundancy(int i);
	void Anneal_Schedule();
	float Calculate_Anneal_Cost();
	void assignLocation(BMP *image);
	void assignLocationType1(list< node>::iterator it, int scale, BMP *image, RGBApixel color);
	void assignLocationType2(int opertype, list< node>::iterator it, list< node>::iterator it_tmp, int scale, BMP *image, RGBApixel color);
	void assignLocationType3(int opertype, list< node>::iterator it, list< node>::iterator it_tmp, int scale, BMP *image, RGBApixel color);
	void assignLocationType4(int opertype, list< node>::iterator it, list< node>::iterator it_tmp, int scale, BMP *image, RGBApixel color);
	void Anneal_move1();
	void Anneal_move2();
	void Anneal_move3();
	void Polish_Exp_Update(bool ctrl);
	void FloorPLot(int *NUM);
	void plot_module(BMP *image);	
	void printPolish_Exp();
	void write_table();
	
};
#endif