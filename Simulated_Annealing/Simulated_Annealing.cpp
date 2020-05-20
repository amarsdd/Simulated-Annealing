#include<iostream>
#include<random>
#include<vector>
#include<list>
#include<iterator>
#include<fstream>
#include<string>
#include"Simulated_Annealing.h"
#include"EasyBMP.h"
#include"EasyBMP_Geometry.h"
#include "EasyBMP_Font.h"

using namespace std;
uniform_real_distribution<> dis(0, 1);
mt19937 Rand_gen;

void initialize_rand()
{
	Rand_gen.seed(time(NULL));
}

int get_Random(int min, int max)
{
	float rnd;
	rnd = dis(Rand_gen)*(max - min + 1) + min;
	return (int)rnd;
}

void node::initialize()
{
	this->FinalBdnCurveIdx1 = 0;
	this->FinalBdnCurveIdx2 = 0;
	this->final_node1 = 0;
	this->final_node2 = 0;
	this->Module1 = NULL;
	this->Module2 = NULL;
	this->Operator1 = NULL;
	this->Operator2 = NULL;
	this->next = NULL;
	this->type_flag = 0;
	this->head_flag = 0;
	this->operator_type = 0;
}

void Sim_Annealing::combineBdnCurvesType1(node* temp, int seq_id1, int seq_id2, int i, int operatorflag)
{
	temp->node_type = 1;
	temp->Module2 = &Polish_Exp.ModuleSeq.at(seq_id2 - i - 1);
	temp->Module1 = &Polish_Exp.ModuleSeq.at(seq_id1 - i - 1);
	for (int j = 0; j < temp->Module1->BdnCurveIdx.size(); j++)
	{
		temp->FinalBdnCurveIdx1 = j;
		for (int k = 0; k < temp->Module2->BdnCurveIdx.size(); k++)
		{
			temp->FinalBdnCurveIdx2 = k;
			if (operatorflag == -1)
			{
				temp->w = temp->Module1->BdnCurveIdx.at(j).w + temp->Module2->BdnCurveIdx.at(k).w;
				temp->h = max(temp->Module1->BdnCurveIdx.at(j).h, temp->Module2->BdnCurveIdx.at(k).h);
			}
			else
			{
				temp->w = max(temp->Module1->BdnCurveIdx.at(j).w, temp->Module2->BdnCurveIdx.at(k).w);
				temp->h = temp->Module1->BdnCurveIdx.at(j).h + temp->Module2->BdnCurveIdx.at(k).h;
			}
			Polish_Exp.OperatorSeq.at(i).combineBdnCurve.push_back(*temp);
		}
	}
}

void Sim_Annealing::combineBdnCurvesType2(node* temp, int seq_id1, int seq_id2, int i, int operatorflag)
{
	temp->node_type = 2;
	int a;
	for (a = 0; a < (Polish_Exp.N - 1); a++)
	{
		if (Polish_Exp.ModuleSeq.at(a).Original_Module_ID == Polish_Exp.Full_Polish_Exp.at(seq_id1))
			break;
	}
	temp->Module1 = &Polish_Exp.ModuleSeq.at(a);
	temp->Operator1 = &Polish_Exp.OperatorSeq.at(i - 1);
	list< node>::iterator temp_it;
	for (int j = 0; j < temp->Module1->BdnCurveIdx.size(); j++)
	{
		temp_it = temp->Operator1->combineBdnCurve.begin();
		temp->FinalBdnCurveIdx1 = j;
		for (int k = 0; k < temp->Operator1->combineBdnCurve.size(); k++)
		{
			if (operatorflag == -1)
			{
				temp->w = temp->Module1->BdnCurveIdx.at(j).w + temp_it->w;
				temp->h = max(temp->Module1->BdnCurveIdx.at(j).h, temp_it->h);
			}
			else
			{
				temp->w = max(temp->Module1->BdnCurveIdx.at(j).w, temp_it->w);
				temp->h = temp->Module1->BdnCurveIdx.at(j).h + temp_it->h;
			}
			temp->final_node1 = k;
			Polish_Exp.OperatorSeq.at(i).combineBdnCurve.push_back(*temp);
			temp_it++;
		}
	}
}

void Sim_Annealing::combineBdnCurvesType3(node* temp, int seq_id1, int seq_id2, int i, int operatorflag)
{
	temp->node_type = 3;
	int a;
	for (a = 0; a < (Polish_Exp.N - 2); a++)
	{
		if (Polish_Exp.OperatorSeq.at(a).inPolish_ID == (seq_id1 + 1))
			break;
	}
	temp->Operator1 = &Polish_Exp.OperatorSeq.at(a);
	temp->Operator2 = &Polish_Exp.OperatorSeq.at(i - 1);
	list< node>::iterator it1 = temp->Operator1->combineBdnCurve.begin();
	list< node>::iterator it2;
	for (int j = 0; j < temp->Operator1->combineBdnCurve.size(); j++)
	{
		it2 = temp->Operator2->combineBdnCurve.begin();
		temp->final_node1 = j;
		for (int k = 0; k < temp->Operator2->combineBdnCurve.size(); k++)
		{
			if (operatorflag == -1)
			{
				temp->w = it1->w + it2->w;
				temp->h = max(it1->h, it2->h);
			}
			else
			{
				temp->w = max(it1->w, it2->w);
				temp->h = it1->h + it2->h;
			}
			temp->final_node2 = k;
			Polish_Exp.OperatorSeq.at(i).combineBdnCurve.push_back(*temp);
			it2++;
		}
		it1++;
	}
}
void Sim_Annealing::combineBdnCurvesType4(node* temp, int seq_id1, int seq_id2, int i, int operatorflag)
{
	temp->node_type = 4;
	temp->Module1 = &Polish_Exp.ModuleSeq.at(seq_id2 - i - 1);
	temp->Operator1 = &Polish_Exp.OperatorSeq.at(i - 1);
	list< node>::iterator it;
	for (int j = 0; j < temp->Module1->BdnCurveIdx.size(); j++)
	{
		it = temp->Operator1->combineBdnCurve.begin();
		temp->FinalBdnCurveIdx1 = j;
		for (int k = 0; k < temp->Operator1->combineBdnCurve.size(); k++)
		{
			if (operatorflag == -1)
			{
				temp->w = temp->Module1->BdnCurveIdx.at(j).w + it->w;
				temp->h = max(temp->Module1->BdnCurveIdx.at(j).h, it->h);
			}
			else
			{
				temp->w = max(temp->Module1->BdnCurveIdx.at(j).w, it->w);
				temp->h = temp->Module1->BdnCurveIdx.at(j).h + it->h;
			}
			temp->final_node1 = k;
			Polish_Exp.OperatorSeq.at(i).combineBdnCurve.push_back(*temp);
			it++;
		}
	}
}



float Sim_Annealing::Calculate_Anneal_Cost()
{

	int seq_id2, seq_id1;
	int i = 0;
	node temp;
	temp.node_type = 0;
	for (i; i < (Polish_Exp.N - 1); i++)
	{
		Polish_Exp.OperatorSeq.at(i).combineBdnCurve.clear();
		temp.initialize();
		seq_id2 = Polish_Exp.OperatorSeq.at(i).inPolish_ID - 1;
		int operatorflag = Polish_Exp.OperatorSeq.at(i).Type;

		temp.operator_type = operatorflag;
		if (Polish_Exp.Full_Polish_Exp.at(seq_id2 - 1) > 0)
		{
			seq_id1 = seq_id2 - 1;
			if (Polish_Exp.Full_Polish_Exp.at(seq_id1 - 1) > 0)
			{
				combineBdnCurvesType1(&temp, seq_id1, seq_id2, i, operatorflag);
			}
			else
			{
				combineBdnCurvesType4(&temp, seq_id1, seq_id2, i, operatorflag);
			}
		}
		else
		{
			list< node>::iterator it;
			temp.node_type = 2;
			temp.Operator1 = &Polish_Exp.OperatorSeq.at(i - 1);
			it = temp.Operator1->combineBdnCurve.begin();
			int seq_delt = 0;
			if (Polish_Exp.Full_Polish_Exp.at(seq_id1) > 0)
			{
				combineBdnCurvesType2(&temp, seq_id1, seq_id2, i, operatorflag);
			}
			else
			{
				combineBdnCurvesType3(&temp, seq_id1, seq_id2, i, operatorflag);
			}
		}
		Remove_Redundancy(i);
	}

	return getAreaCost();
}

float Sim_Annealing::getAreaCost()
{
	list< node>::iterator temp_node, smallest_area, smallest_valid_area_pq, smallest_valid_area_pq1;
	temp_node = smallest_valid_area_pq = smallest_valid_area_pq1 = smallest_area = Polish_Exp.OperatorSeq.back().combineBdnCurve.begin();
	float aspect_ratio = temp_node->h / temp_node->w;
	float area = temp_node->h*temp_node->w;
	float valid_area_pq = area;
	float valid_area_pq1 = area;
	for (int x = 1; x < Polish_Exp.OperatorSeq.back().combineBdnCurve.size(); x++)
	{
		temp_node++;
		float area_comp = temp_node->h*temp_node->w;
		aspect_ratio = temp_node->h / temp_node->w;
		/*if ((aspect_ratio >= Polish_Exp.p) && (aspect_ratio <= Polish_Exp.q))
		{
		smallest_valid_area_pq1 = temp_node;
		valid_area_pq1 = area_comp;
		}*/
		if (area_comp < area)
		{
			smallest_area = temp_node;
			area = area_comp;
			if ((aspect_ratio >= Polish_Exp.p) && (aspect_ratio <= Polish_Exp.q))
			{
				smallest_valid_area_pq = temp_node;
				valid_area_pq = area_comp;
			}
		}
	}
	aspect_ratio = smallest_valid_area_pq->h / smallest_valid_area_pq->w;
	Polish_Exp.final_node = smallest_area;

	if ((aspect_ratio >= Polish_Exp.p) && (aspect_ratio <= Polish_Exp.q))
	{
		return valid_area_pq;
	}
	else
	{
		return 2 * (area);
	}
}
void Sim_Annealing::Remove_Redundancy(int i)
{
	list< node>::iterator it1, it2;
	float h1, w1, h2, w2;
	it1 = Polish_Exp.OperatorSeq.at(i).combineBdnCurve.begin();
	while (it1 != Polish_Exp.OperatorSeq.at(i).combineBdnCurve.end())
	{
		int check = 0;
		h1 = it1->h; 
		w1 = it1->w;
		it2 = it1;
		it2++;
		while (it2 != Polish_Exp.OperatorSeq.at(i).combineBdnCurve.end())
		{
			h2 = it2->h; 
			w2 = it2->w;
			if ((h1 >= h2) && (w1 >= w2))
			{
				it1 = Polish_Exp.OperatorSeq.at(i).combineBdnCurve.erase(it1);
				check = 1;
				break;
			}
			else
			{
				if ((h1 <= h2) && (w1 <= w2))
				{
					it2 = Polish_Exp.OperatorSeq.at(i).combineBdnCurve.erase(it2);
					continue;
				}
			}
			it2++;
		}
	}
}

void Sim_Annealing::Anneal_move1()
{
	int rand = get_Random(0, Polish_Exp.N - 2);
	int id1 = Polish_Exp.ModuleSeq.at(rand).inPolish_ID;
	int id2 = Polish_Exp.ModuleSeq.at(rand + 1).inPolish_ID;
	Polish_Exp.ModuleSeq.at(rand).inPolish_ID = id2;
	Polish_Exp.ModuleSeq.at(rand + 1).inPolish_ID = id1;
	vector<Module>::iterator first = Polish_Exp.ModuleSeq.begin();
	vector<Module>::iterator second = Polish_Exp.ModuleSeq.begin();
	iter_swap(first + rand, second + rand + 1);
	Polish_Exp_Update(true);
}
void Sim_Annealing::Anneal_move2()
{
	int rand = get_Random(0, Polish_Exp.Operator_chain.size() - 1);
	list<Operator_chain>::iterator it = Polish_Exp.Operator_chain.begin();
	advance(it, rand);
	int p = it->Operator_Serial;
	int pos = Polish_Exp.OperatorSeq.at(p).inPolish_ID - 1;
	while (pos < (2 * Polish_Exp.N - 1))
	{
		if (Polish_Exp.Full_Polish_Exp.at(pos) == -1)
		{
			Polish_Exp.OperatorSeq.at(p).Type = -2;
			p++;
		}
		else
		{
			if (Polish_Exp.Full_Polish_Exp.at(pos) == -2)
			{
				Polish_Exp.OperatorSeq.at(p).Type = -1;
				p++;
			}
			else
				break;
		}
		pos++;
	}
	Polish_Exp_Update(true);
}
void Sim_Annealing::Anneal_move3()
{
	bool flag = false;
	int count = 0;
	while ((flag != true) && count < (10 * Polish_Exp.N))
	{
		int randchoice = get_Random(1, 2);
		int randpos;
		list<Operator_chain>::iterator it = Polish_Exp.Operator_chain.begin();
		switch (randchoice)
		{
		case 1:
			randpos = get_Random(0, Polish_Exp.Operator_chain.size() - 1);
			advance(it, randpos);
			if (Polish_Exp.Full_Polish_Exp.at(it->inPolish_ID - 2) != Polish_Exp.Full_Polish_Exp.at(it->inPolish_ID))
			{
				int dk = 2 * (it->Operator_Serial + 1);
				if (dk < it->inPolish_ID)
				{
					int id1, id2;
					id1 = Polish_Exp.ModuleSeq.at(it->Module_Serial).inPolish_ID;
					id2 = Polish_Exp.OperatorSeq.at(it->Operator_Serial).inPolish_ID;
					Polish_Exp.ModuleSeq.at(it->Module_Serial).inPolish_ID = id2;
					Polish_Exp.OperatorSeq.at(it->Operator_Serial).inPolish_ID = id1;
					flag = true;
				}
			}
			break;
		case 2:
			if (Polish_Exp.Operator_chain.size() < 2)
				break;
			randpos = get_Random(0, Polish_Exp.Operator_chain.size() - 2);
			advance(it, randpos);
			if (Polish_Exp.Full_Polish_Exp.at(it->inPolish_ID_R) != Polish_Exp.Full_Polish_Exp.at(it->inPolish_ID_R + 2))
			{
				int id1, id2;
				id1 = Polish_Exp.ModuleSeq.at(it->Module_Serial_R).inPolish_ID;
				id2 = Polish_Exp.OperatorSeq.at(it->Operator_Serial_R).inPolish_ID;
				Polish_Exp.ModuleSeq.at(it->Module_Serial_R).inPolish_ID = id2;
				Polish_Exp.OperatorSeq.at(it->Operator_Serial_R).inPolish_ID = id1;
				flag = true;
			}
			break;
		}
		count++;
	}
	if (flag == true)
	{
		Polish_Exp_Update(true);
	}
}


void Sim_Annealing::Anneal_Schedule()
{
	float accept_count;
	float accept_ratio = 0;
	float cost, delta_cost;
	FloorplanData data_table;
	Polish_Expression Polish_Expressionbackup = Polish_Exp;
	while ((T > T_final) && (accept_ratio > 0.01))
	{
		int i;
		accept_count = 0;
		data_table.area = Polish_Exp.final_node->h*Polish_Exp.final_node->w;
		data_table.t = T;
		table.push_back(data_table);
		for (i = 0; i < numMoves; i++)
		{
			int rand = get_Random(1, 3);
			switch (rand)
			{
			case 1:
				Anneal_move1();
				num_move1++;
				break;
			case 2:
				Anneal_move2();
				num_move2++;
				break;
			case 3:
				Anneal_move3();
				num_move3++;
				break;
			}
			cost = Calculate_Anneal_Cost();
			delta_cost = cost - last_accepted_cost;
			float rand01 = dis(Rand_gen);
			if ((delta_cost < 0) || (rand01 < exp(-delta_cost / T)))
			{
				Polish_Expressionbackup = Polish_Exp;
				last_accepted_cost = cost;
				accept_count++;
			}
			else
			{
				Polish_Exp = Polish_Expressionbackup;
			}
			if (delta_cost < 0)
			{
				Polish_Expressionbest = Polish_Exp;
			}
		}
		Calculate_Anneal_Cost();
		accept_ratio = accept_count / numMoves;
		T = T_dec_ratio*T;
	}
	Polish_Exp = Polish_Expressionbest;
	Calculate_Anneal_Cost();
}


void Sim_Annealing::Init_Polish_Exp()
{
	Polish_Exp.ModuleSeq.at(0).inPolish_ID = 1;
	Polish_Exp.ModuleSeq.at(0).inPolish_ID = 1;
	Polish_Exp.ModuleSeq.at(0).FinalBdnCurveIdx = 0;
	Polish_Exp.ModuleSeq.at(0).FinalBdnCurveIdx_tmp = 0;
	int seq, i, operatype;
	float net_area_val = 0;
	seq = 2;
	operatype = -1;
	net_area_val += Polish_Exp.ModuleSeq.at(0).A_;
	BdnCurveIdxCal(&Polish_Exp.ModuleSeq.at(0));
	for (i = 1; i < Polish_Exp.N; i++)
	{
		Polish_Exp.ModuleSeq.at(i).inPolish_ID = seq++;
		Polish_Exp.ModuleSeq.at(i).FinalBdnCurveIdx = 0;
		Polish_Exp.ModuleSeq.at(i).FinalBdnCurveIdx_tmp = 0;
		net_area_val += Polish_Exp.ModuleSeq.at(i).A_;
		BdnCurveIdxCal(&Polish_Exp.ModuleSeq.at(i));
		Polish_Exp.OperatorSeq.at(i - 1).inPolish_ID = seq++;
		Polish_Exp.OperatorSeq.at(i - 1).Type = operatype;
		if (operatype == -1)
			operatype = -2;
		else
			operatype = -1;
	}
	Polish_Exp_Update(true);
}
void Sim_Annealing::Polish_Exp_Update(bool ctrl)
{

	for (int i = 0; i < Polish_Exp.N; i++)
	{
		Polish_Exp.Full_Polish_Exp.at(Polish_Exp.ModuleSeq.at(i).inPolish_ID - 1) = Polish_Exp.ModuleSeq.at(i).Original_Module_ID;
		if (i == (Polish_Exp.N - 1))
			continue;
		Polish_Exp.Full_Polish_Exp.at(Polish_Exp.OperatorSeq.at(i).inPolish_ID - 1) = Polish_Exp.OperatorSeq.at(i).Type;
	}
	if (ctrl == false)
		return;
	Polish_Exp.Operator_chain.clear();
	Operator_chain Operatorchain;
	vector<int>::iterator it = Polish_Exp.Full_Polish_Exp.begin();
	int count = 0, operatorcount = 0, Modulecount = 0;
	while (it != Polish_Exp.Full_Polish_Exp.end())
	{
		while (*it > 0)
		{
			it++;
			count++;
			Modulecount++;
		}
		Operatorchain.inPolish_ID = count;
		Operatorchain.Module_Serial = Modulecount - 1;
		Operatorchain.Operator_Serial = operatorcount;
		Operatorchain.inPolish_ID_R = -1;
		Operatorchain.Module_Serial_R = -1;
		Operatorchain.Operator_Serial_R = -1;
		Polish_Exp.Operator_chain.push_back(Operatorchain);
		while (*it < 0)
		{
			it++;
			count++;
			operatorcount++;
			if (it == Polish_Exp.Full_Polish_Exp.end())
				break;
		}
		if (it != Polish_Exp.Full_Polish_Exp.end())
		{
			Polish_Exp.Operator_chain.back().inPolish_ID_R = count - 1;
			Polish_Exp.Operator_chain.back().Module_Serial_R = Modulecount;
			Polish_Exp.Operator_chain.back().Operator_Serial_R = operatorcount - 1;
		}
	}
}
void Sim_Annealing::BdnCurveIdxCal(Module *Module)
{
	int i = 0;
	float delta, count = 0;
	BdnCurve BdnCurveIdx;
	delta = sqrt(Module->A_ / Module->r_) - sqrt(Module->A_ / Module->s_);
	delta /= (interpol_num - 1);
	while (i < interpol_num)
	{
		BdnCurveIdx.w = count*delta + sqrt(Module->A_ / Module->s_);
		BdnCurveIdx.h = Module->A_ / BdnCurveIdx.w;
		BdnCurveIdx.BdnCurveID = i;
		Module->BdnCurveIdx.push_back(BdnCurveIdx);
		count++;
		i++;
		if (Module->r_ == Module->s_)
		{
			break;
		}
	}
	if (Module->S == 1)
		return;
	delta = sqrt(Module->A_ * Module->s_) - sqrt(Module->A_ * Module->r_);
	delta /= (interpol_num - 1);
	count = 0;
	while (i < 2 * interpol_num)
	{
		BdnCurveIdx.w = count*delta + sqrt(Module->A_ * Module->r_);
		BdnCurveIdx.h = Module->A_ / BdnCurveIdx.w;
		BdnCurveIdx.BdnCurveID = i;
		Module->BdnCurveIdx.push_back(BdnCurveIdx);
		count++;
		i++;
		if (Module->r_ == Module->s_)
		{
			break;
		}
	}
}

void Sim_Annealing::assignLocationType1(list< node>::iterator it, int scale, BMP *image, RGBApixel color)
{
	int opertype;
	opertype = it->operator_type;
	it->Module1->FinalBdnCurveIdx = it->FinalBdnCurveIdx1;
	it->Module2->FinalBdnCurveIdx = it->FinalBdnCurveIdx2;

	it->Module1->parentnode = &(*it);
	it->Module2->parentnode = &(*it);

	it->Module2->location[0] = it->location[0];
	it->Module2->location[1] = it->location[1];
	if (opertype == -1)
	{
		it->Module1->location[0] = it->location[0] - it->Module2->BdnCurveIdx.at(it->FinalBdnCurveIdx2).w;
		it->Module1->location[1] = it->location[1];
		it->Module2->center[0] = it->location[0] - it->Module2->BdnCurveIdx.at(it->FinalBdnCurveIdx2).w / 2;
		it->Module2->center[1] = it->location[1] - it->h / 2;
		it->Module1->center[0] = it->Module1->location[0] - it->Module1->BdnCurveIdx.at(it->FinalBdnCurveIdx1).w / 2;
		it->Module1->center[1] = it->location[1] - it->h / 2;
	}
	else
	{
		it->Module1->location[1] = it->location[1] - it->Module2->BdnCurveIdx.at(it->FinalBdnCurveIdx2).h;
		it->Module1->location[0] = it->location[0];
		it->Module2->center[1] = it->location[1] - it->Module2->BdnCurveIdx.at(it->FinalBdnCurveIdx2).h / 2;
		it->Module2->center[0] = it->location[0] - it->w / 2;
		it->Module1->center[1] = it->Module1->location[1] - it->Module1->BdnCurveIdx.at(it->FinalBdnCurveIdx1).h / 2;
		it->Module1->center[0] = it->location[0] - it->w / 2;
	}
	it->type_flag = 3;

}

void Sim_Annealing::assignLocationType2(int opertype, list< node>::iterator it, list< node>::iterator it_tmp, int scale, BMP *image, RGBApixel color)
{

	it->location[0] = it_tmp->location[0];
	it->location[1] = it_tmp->location[1];
	if (opertype == -1)
	{
		it_tmp->Module1->location[0] = it_tmp->location[0] - it->w;
		it_tmp->Module1->location[1] = it_tmp->location[1];
		it_tmp->Module1->center[0] = it_tmp->Module1->location[0] - it_tmp->Module1->BdnCurveIdx.at(it_tmp->FinalBdnCurveIdx1).w / 2;
		it_tmp->Module1->center[1] = it_tmp->location[1] - it_tmp->h / 2;
	}
	else
	{
		it_tmp->Module1->location[1] = it_tmp->location[1] - it->h;
		it_tmp->Module1->location[0] = it_tmp->location[0];
		it_tmp->Module1->center[1] = it_tmp->Module1->location[1] - it_tmp->Module1->BdnCurveIdx.at(it_tmp->FinalBdnCurveIdx1).h / 2;
		it_tmp->Module1->center[0] = it_tmp->location[0] - it_tmp->w / 2;
	}
	it->parentnode = it_tmp;
}

void Sim_Annealing::assignLocationType3(int opertype, list< node>::iterator it, list< node>::iterator it_tmp, int scale, BMP *image, RGBApixel color)
{
	if (it->type_flag == 0)
	{
		list< node>::iterator it_temp;
		it->type_flag = 1;
		opertype = it->operator_type;
		it_tmp = it;
		it_temp = it_tmp->Operator1->combineBdnCurve.begin();
		advance(it_temp, it_tmp->final_node1);
		it = it->Operator2->combineBdnCurve.begin();
		advance(it, it_tmp->final_node2);
		it->location[0] = it_tmp->location[0];
		it->location[1] = it_tmp->location[1];
		it->parentnode = it_tmp;
	}
	else
	{
		opertype = it->operator_type;
		list< node>::iterator it_temp;
		it->type_flag = 3;
		it_tmp = it;
		it = it->Operator1->combineBdnCurve.begin();
		advance(it, it_tmp->final_node1);
		it_temp = it_tmp->Operator2->combineBdnCurve.begin();
		advance(it_temp, it_tmp->final_node2);
		if (opertype == -1)
		{
			it->location[0] = it_tmp->location[0] - it_temp->w;
			it->location[1] = it_tmp->location[1];

		}
		else
		{
			it->location[1] = it_tmp->location[1] - it_temp->h;
			it->location[0] = it_tmp->location[0];

		}
		it->parentnode = it_tmp;
	}
}

void Sim_Annealing::assignLocationType4(int opertype, list< node>::iterator it, list< node>::iterator it_tmp, int scale, BMP *image, RGBApixel color)
{
	it_tmp->Module1->location[0] = it_tmp->location[0];
	it_tmp->Module1->location[1] = it_tmp->location[1];
	if (opertype == -1)
	{
		it->location[0] = it_tmp->location[0] - it_tmp->Module1->BdnCurveIdx.at(it_tmp->FinalBdnCurveIdx1).w;
		it->location[1] = it_tmp->location[1];
		it_tmp->Module1->center[0] = it_tmp->location[0] - it_tmp->Module1->BdnCurveIdx.at(it_tmp->FinalBdnCurveIdx1).w / 2;
		it_tmp->Module1->center[1] = it_tmp->location[1] - it_tmp->h / 2;
	}
	else
	{
		it->location[0] = it_tmp->location[0];
		it->location[1] = it_tmp->location[1] - it_tmp->Module1->BdnCurveIdx.at(it_tmp->FinalBdnCurveIdx1).h;
		it_tmp->Module1->center[1] = it_tmp->location[1] - it_tmp->Module1->BdnCurveIdx.at(it_tmp->FinalBdnCurveIdx1).h / 2;
		it_tmp->Module1->center[0] = it_tmp->location[0] - it_tmp->w / 2;
	}
	it->parentnode = it_tmp;
}

void Sim_Annealing::assignLocation(BMP *image)
{
	float scale = 50;
	float x1, x2, y2, y1;
	RGBApixel color;
	color.Blue = 0;
	color.Green = 0;
	color.Red = 0;
	list< node>::iterator it, it_tmp;
	Polish_Exp.final_node->head_flag = 1;
	Polish_Exp.final_node->location[0] = Polish_Exp.final_node->w;
	Polish_Exp.final_node->location[1] = Polish_Exp.final_node->h;
	int opertype = Polish_Exp.final_node->operator_type;
	it_tmp = it = Polish_Exp.final_node;
	int nodesel, leafsel;
	while (1)
	{
		switch (it->node_type)
		{
		case 1:
			assignLocationType1(it, scale, image, color);
			if (it->head_flag == 1)
				return;
			it = it->parentnode;
			break;
		case 2:
			if (it->type_flag == 3)
			{
				if (it->head_flag != 1)
					it = it->parentnode;
				else
					return;
			}
			else
			{
				opertype = it->operator_type;
				it->Module1->FinalBdnCurveIdx = it->FinalBdnCurveIdx1;
				it->Module1->parentnode = &(*it);
				it->type_flag = 3;
				it_tmp = it;
				it = it->Operator1->combineBdnCurve.begin();
				advance(it, it_tmp->final_node1);
				assignLocationType2(opertype, it, it_tmp, scale, image, color);
			}
			break;
		case 3:
			if (it->type_flag == 3)
			{
				if (it->head_flag != 1)
					it = it->parentnode;
				else
					return;
			}
			else
			{
				assignLocationType2(opertype, it, it_tmp, scale, image, color);
			}
			break;
		case 4:
			if (it->type_flag == 3)
			{
				if (it->head_flag != 1)
					it = it->parentnode;
				else
					return;
			}
			else
			{
				opertype = it->operator_type;
				it->type_flag = 3;
				nodesel = it->final_node1;
				leafsel = it->FinalBdnCurveIdx1;
				it->Module1->FinalBdnCurveIdx = leafsel;
				it->Module1->parentnode = &(*it);
				it_tmp = it;
				it = it->Operator1->combineBdnCurve.begin();
				advance(it, nodesel);
				assignLocationType4(opertype, it, it_tmp, scale, image, color);
			}
			break;
		}
	}
}

void Sim_Annealing::printPolish_Exp()
{
	for (int i = 0; i < Polish_Exp.Full_Polish_Exp.size(); i++)
	{
		if (Polish_Exp.Full_Polish_Exp.at(i)>0)
			printf("%d ", Polish_Exp.Full_Polish_Exp.at(i));
		else
		{
			if (Polish_Exp.Full_Polish_Exp.at(i) == -1)
				cout << "* ";
			else
				cout << "+ ";

		}
	}
	cout << "\n";
}
void Sim_Annealing::FloorPLot(int *NUM)
{
	BMP image;
	RGBApixel color;
	float scale = 50;
	color.Blue = 0;
	color.Green = 0;
	color.Red = 0;
	image.SetSize(Polish_Exp.final_node->w * scale, Polish_Exp.final_node->h * scale);
	DrawFastLine(image, 0, 0, 0, Polish_Exp.final_node->h * scale, color);
	DrawFastLine(image, 0, 0, Polish_Exp.final_node->w * scale, 0, color);
	DrawFastLine(image, Polish_Exp.final_node->w * scale, Polish_Exp.final_node->h * scale, 0, Polish_Exp.final_node->h * scale, color);
	DrawFastLine(image, Polish_Exp.final_node->w * scale, Polish_Exp.final_node->h * scale, Polish_Exp.final_node->w * scale, 0, color);
	assignLocation(&image);
	plot_module(&image);
	char filename[32];
	sprintf_s(filename, "FloorPlot_%d.bmp", *NUM);
	image.WriteToFile(filename);
}

void Sim_Annealing::plot_module(BMP *image)
{
	char buf[32];
	RGBApixel color;
	float scale = 50;
	color.Blue = 0;
	color.Green = 0;
	color.Red = 0;
	RGBApixel color2, color3;
	color2.Blue = 0;
	color2.Green = 0;
	color2.Red = 255;
	color3.Blue = 255;
	color3.Green = 0;
	color3.Red = 50;
	float x0, x1, y0, y1;
	int finalsel;
	for (int i = 0; i < Polish_Exp.N; i++)
	{
		finalsel = Polish_Exp.ModuleSeq.at(i).FinalBdnCurveIdx;
		y0 = (scale * (Polish_Exp.ModuleSeq.at(i).location[1] - Polish_Exp.ModuleSeq.at(i).BdnCurveIdx[finalsel].h));
		y1 = (scale * (Polish_Exp.ModuleSeq.at(i).location[1]));
		x0 = (scale * (Polish_Exp.ModuleSeq.at(i).location[0] - Polish_Exp.ModuleSeq.at(i).BdnCurveIdx[finalsel].w));
		x1 = (scale * (Polish_Exp.ModuleSeq.at(i).location[0]));

		/*float diff_h = (Polish_Exp.ModuleSeq.at(i).location[1] - Polish_Exp.ModuleSeq.at(i).center[1]) * 2 - Polish_Exp.ModuleSeq.at(i).BdnCurveIdx[finalsel].h;
		float diff_w = (Polish_Exp.ModuleSeq.at(i).location[0] - Polish_Exp.ModuleSeq.at(i).center[0]) * 2 - Polish_Exp.ModuleSeq.at(i).BdnCurveIdx[finalsel].w;
		*/
		/*float diff_h = (Polish_Exp.ModuleSeq.at(i).parentnode->h) - Polish_Exp.ModuleSeq.at(i).BdnCurveIdx[finalsel].h;
		float diff_w = (Polish_Exp.ModuleSeq.at(i).parentnode->w) - Polish_Exp.ModuleSeq.at(i).BdnCurveIdx[finalsel].w;

		float xe0, xe1, ye0, ye1;
		if (diff_h > 0 && Polish_Exp.ModuleSeq.at(i).parentnode->operator_type == (-1))
		{
		ye0 = (scale * (Polish_Exp.ModuleSeq.at(i).location[1] - Polish_Exp.ModuleSeq.at(i).BdnCurveIdx[finalsel].h - diff_h));
		ye1 = (scale * (Polish_Exp.ModuleSeq.at(i).location[1] - Polish_Exp.ModuleSeq.at(i).BdnCurveIdx[finalsel].h));
		xe0 = (scale * (Polish_Exp.ModuleSeq.at(i).location[0] - Polish_Exp.ModuleSeq.at(i).BdnCurveIdx[finalsel].w));
		xe1 = (scale * (Polish_Exp.ModuleSeq.at(i).location[0]));
		DrawFastLine(*image, xe0, ye0, xe1, ye0, color);
		DrawFastLine(*image, xe0, ye0, xe0, ye1, color);
		DrawFastLine(*image, xe1, ye1, xe0, ye1, color);
		DrawFastLine(*image, xe1, ye1, xe1, ye0, color);
		DrawFastLine(*image, xe1, ye1, xe0, ye0, color2);
		DrawFastLine(*image, xe1, ye0, xe0, ye1, color2);
		}
		else if (diff_w > 0 && Polish_Exp.ModuleSeq.at(i).parentnode->operator_type == (-2))
		{
		ye0 = (scale * (Polish_Exp.ModuleSeq.at(i).location[1] - Polish_Exp.ModuleSeq.at(i).BdnCurveIdx[finalsel].h));
		ye1 = (scale * (Polish_Exp.ModuleSeq.at(i).location[1]));
		xe0 = (scale * (Polish_Exp.ModuleSeq.at(i).location[0] - Polish_Exp.ModuleSeq.at(i).BdnCurveIdx[finalsel].w - diff_w));
		xe1 = (scale * (Polish_Exp.ModuleSeq.at(i).location[0] - Polish_Exp.ModuleSeq.at(i).BdnCurveIdx[finalsel].w));
		DrawFastLine(*image, xe0, ye0, xe1, ye0, color);
		DrawFastLine(*image, xe0, ye0, xe0, ye1, color);
		DrawFastLine(*image, xe1, ye1, xe0, ye1, color);
		DrawFastLine(*image, xe1, ye1, xe1, ye0, color);
		DrawFastLine(*image, xe1, ye1, xe0, ye0, color2);
		DrawFastLine(*image, xe1, ye0, xe0, ye1, color2);
		}*/

		DrawFastLine(*image, x0, y0, x1, y0, color);
		DrawFastLine(*image, x0, y0, x0, y1, color);
		DrawFastLine(*image, x1, y1, x0, y1, color);
		DrawFastLine(*image, x1, y1, x1, y0, color);

		if (Polish_Exp.ModuleSeq.at(i).parentnode->operator_type == (-1))
		{
			int div = (int)((x1 - x0) / 10);

			for (float xt = x0 + div; xt < x1; xt = xt + div)
			{
				DrawFastLine(*image, xt, y0, xt, y1, color2);
			}
			for (float yt = y0 + div; yt < y1; yt = yt + div)
			{
				DrawFastLine(*image, x0, yt, x1, yt, color2);
			}
		}
		if (Polish_Exp.ModuleSeq.at(i).parentnode->operator_type == (-2))
		{
			int div = (int)((y1 - y0) / 10);

			for (float xt = x0 + div; xt < x1; xt = xt + div)
			{
				DrawFastLine(*image, xt, y0, xt, y1, color2);
			}
			for (float yt = y0 + div; yt < y1; yt = yt + div)
			{
				DrawFastLine(*image, x0, yt, x1, yt, color2);
			}
		}
		sprintf_s(buf, "M%d\n", Polish_Exp.ModuleSeq.at(i).Original_Module_ID);
		PrintString(*image, buf, (x1 + x0) / 2 - 30, (y1 + y0) / 2 - 20, 30, color3);
	}
}

bool Sim_Annealing::InputfileRead()
{
	tmp_ = 0;
	string filename;
	filename = "HW4.txt";
	ifstream netlist_file(filename);
	if (!netlist_file.is_open())
	{
		printf("file open failed!\n");
		return false;
	}
	string line;
	getline(netlist_file, line);
	string::size_type pos = 0;
	Polish_Exp.N = (int)get_data(line, pos);
	pos = line.find(" ", pos + 1);
	if (pos == string::npos)
	{
		printf("header error!\n");
		return false;
	}
	Polish_Exp.p = get_data(line, pos + 1);
	pos = line.find(" ", pos + 1);
	if (pos == string::npos)
	{
		printf("header error!\n");
		return false;
	}
	Polish_Exp.q = get_data(line, pos + 1);


	Polish_Exp.ModuleSeq.resize(Polish_Exp.N);
	Polish_Exp.OperatorSeq.resize(Polish_Exp.N - 1);
	Polish_Exp.Full_Polish_Exp.resize(2 * Polish_Exp.N - 1, 0);


	while (getline(netlist_file, line))
	{
		if (line.size() < 5)
			continue;
		string::size_type pos = 0;
		int i = 0;
		Polish_Exp.ModuleSeq.at(tmp_).Original_Module_ID = (int)get_data(line, pos);
		pos = line.find(" ", pos + 1);
		if (pos == string::npos)
		{
			printf("body error!\n");

		}
		Polish_Exp.ModuleSeq.at(tmp_).A_ = get_data(line, pos + 1);
		pos = line.find(" ", pos + 1);
		if (pos == string::npos)
		{
			printf("body error!\n");

		}
		Polish_Exp.ModuleSeq.at(tmp_).r_ = get_data(line, pos + 1);
		pos = line.find(" ", pos + 1);
		if (pos == string::npos)
		{
			printf("body error!\n");

		}
		Polish_Exp.ModuleSeq.at(tmp_).s_ = get_data(line, pos + 1);
		pos = line.find(" ", pos + 1);
		if (pos == string::npos)
		{
			printf("body error!\n");

		}
		Polish_Exp.ModuleSeq.at(tmp_).S = get_data(line, pos + 1);

		tmp_++;
	}
	netlist_file.close();
	return true;
}
float Sim_Annealing::get_data(string line, string::size_type pos)
{
	string::size_type tmp;
	string str;
	str.clear();
	for (tmp = pos; tmp < line.size(); tmp++)
	{
		if ((line[tmp] >= '0') && (line[tmp] <= '9') || (line[tmp] == '.'))
			str += line[tmp];
		else
			break;
	}
	return stof(str);
}



void Sim_Annealing::write_table()
{
	ofstream out;
	out.open("tables.txt");
	if (out.is_open())
	{
		for (int i = 0; i < Polish_Exp.N; i++)
		{
			int sel = Polish_Exp.ModuleSeq[i].FinalBdnCurveIdx;
			out << Polish_Exp.ModuleSeq[i].Original_Module_ID << "	" << Polish_Exp.ModuleSeq[i].A_ << "	" << Polish_Exp.ModuleSeq[i].r_ << "	" << Polish_Exp.ModuleSeq[i].s_ << "	" << Polish_Exp.ModuleSeq[i].S << "	" << "	" << Polish_Exp.ModuleSeq[i].BdnCurveIdx[sel].h << "	" << Polish_Exp.ModuleSeq[i].BdnCurveIdx[sel].w << std::endl;

		}
		out << std::endl << std::endl << std::endl;
		for (int i = 0; i < this->table.size(); i++)
		{
			out << this->table[i].t << "	" << this->table[i].area << std::endl;
		}
		out.close();
	}
	else
	{
		cout << "Error opening file to write table";
	}
}


