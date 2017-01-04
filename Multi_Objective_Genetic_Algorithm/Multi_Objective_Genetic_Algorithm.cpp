//

#include "stdafx.h"
#include <iostream>

using std::cout;
using std::endl;

inline float RandomFloat(float min, float max)
{
	float r = (float)rand() / (float)RAND_MAX;
	return min + r * (max - min);
}

inline float RandomInt(int min, int max)
{
	float r = (float)rand() / (float)RAND_MAX;
	return (int)((float)min + r * float(max - min));
}

namespace BIT_OPS
{

	void SetBit(int &dna, int idx)
	{
		dna |= (1 << idx);
	}

	bool CheckBit(int &dna, int idx)
	{
		return dna & (1 << idx);
	}

	void ResetBit(int &dna, int idx)
	{
		dna &= ~(1 << idx);
	}


}

namespace DNA_OUTPUT
{
	using namespace BIT_OPS;
	using std::cout;
	using std::endl;

	void print_dna(int g)
	{
		for (int i = 0; i < 32; i++)
		{
			if (CheckBit(g, i)) std::cout << "1";
			else std::cout << "0";
		}
	}
}



namespace PEAKS_SURFACE
{
	float computePeaks(float X, float Y)
	{
		return 3 * (1 - X)*(1 - X)*std::exp(-(X*X) - (Y + 1)*(Y + 1)) - 10 * (X / 5 - std::pow(X, 3) - std::pow(Y, 5))*std::exp(-X*X - Y*Y) - 1 / 3 * std::exp(-(X + 1)*(X + 1) - Y*Y);
	}

	float peaksDerivX(float x, float y)
	{
		return -6 * (1 - x)*(1 - x)*x*std::exp(-x*x - (y + 1)*(y + 1)) - 6 * (1 - x)*std::exp(-x*x - (y + 1)*(y + 1)) + 20 * x*(-std::pow(x, 3) + x / 5 - std::pow(y, 5))*std::exp(-x*x - y*y) - 10 * (1 / 5 - 3 * x*x)*std::exp(-x*x - y*y) - (-2 * x - 2)*std::exp((-x - 1)*(x + 1) - y*y) / 3;
	}

	float peaksDerivY(float x, float y)
	{
		return -6 * (1 - x)*(1 - x)*(y + 1)*std::exp(-(y + 1)*(y + 1) - x*x) + 20 * y*(-std::pow(y, 5) - std::pow(x, 3) + x / 5)*std::exp(-y*y - x*x) + 50 * std::pow(y, 4)*std::exp(-y*y - x*x) + 2 * y*std::exp((-x - 1)*(x + 1) - y*y) / 3;
	}
}

#include <vector>




namespace SURFACE_CLIMBING_GA_2
{
	using std::vector;

	using namespace BIT_OPS;
	using namespace DNA_OUTPUT;

	struct gene
	{
		int x = 0;
		int y = 0;
	};

	enum GENDER { MALE, FEMALE };
	struct individual
	{
		gene m;
		gene f;

		float fitness;
	};

	vector< individual > males;

	vector< individual > females;

	vector< individual > mating_pool_males;

	vector< individual > mating_pool_females;


	void CreatePopulation(int num_members)
	{
		for (int i = 0; i < num_members; i++)
		{
			individual a;
			males.push_back(a);
			mating_pool_males.push_back(a);

			individual b;
			females.push_back(b);
			mating_pool_females.push_back(b);
		}
	}

	void CopyPopulation(vector<individual> &src, vector<individual> &dst)
	{
		dst.clear();
		for (int i = 0; i < src.size(); i++)
		{
			dst.push_back(src[i]);
		}
	}

	float surface_min_x = -3.0f;
	float surface_min_y = -3.0f;

	float surface_max_x = 3.0f;
	float surface_max_y = 3.0f;


	bool GeneInBounds(gene &g)
	{
		float x = ((float)(g.x)) / 10000;
		float y = ((float)(g.y)) / 10000;

		if ((x >= surface_min_x && x <= surface_max_x) &&
			(y >= surface_min_y && y <= surface_max_y))
		{
			return true;
		}

		return false;
	}

	float GeneHeight_Males(gene& g)
	{
		float x = ((float)(g.x)) / 10000;
		float y = ((float)(g.y)) / 10000;

		return PEAKS_SURFACE::computePeaks(x, y);
	}

	float GeneHeight_Females(gene& g)
	{
		float x = ((float)(g.x)) / 10000;
		float y = ((float)(g.y)) / 10000;

		return PEAKS_SURFACE::computePeaks(y, x);
	}

	int int_crossoverSinglePoint(int dnaA, int dnaB)
	{
		int cross_point = RandomInt(0, 32);

		int output = 0;

		for (int i = 0; i < cross_point; i++)
		{
			if (CheckBit(dnaA, i))
				SetBit(output, i);
		}
		for (int i = cross_point; i < 32; i++)
		{
			if (CheckBit(dnaB, i))
				SetBit(output, i);
		}

		return output;
	}

	void GeneCrossOverSinglePoint(gene& output, gene &a, gene& b)
	{
		output.x = int_crossoverSinglePoint(a.x, b.x);
		output.y = int_crossoverSinglePoint(a.y, b.y);
	}

	float mutation_rate = 0.3f;
	int best_index_males = 0;
	int best_index_females = 0;
	float sumFitness_males = 0.0f;
	float sumFitness_females = 0.0f;


	int int_mutate(int x)
	{
		int X = x;

		float mutation_prob = RandomFloat(0, 1.0f);
		if (mutation_prob < mutation_rate)
		{
			int mut_point = RandomInt(5, 15);

			if (CheckBit(X, mut_point))
			{
				ResetBit(X, mut_point);
			}
			else
			{
				SetBit(X, mut_point);
			}
		}

		return X;
	}

	void GeneMutate(gene& g)
	{
		g.x = int_mutate(g.x);
		g.y = int_mutate(g.y);
	}


	void MutatePopulation()
	{
		for (int i = 0; i < males.size(); i++)
		{
			if (i != best_index_males)
			{
				GeneMutate(males[i].m);
				GeneMutate(males[i].f);
			}
			if (i != best_index_males)
			{
				GeneMutate(females[i].m);
				GeneMutate(females[i].f);
			}
		}
	}

	void RamdomizeInitialPopulation()
	{
		for (int i = 0; i < females.size(); i++)
		{
			males[i].m.x = (int)(RandomFloat(0, 60000.0) - 30000.0);
			males[i].m.y = (int)(RandomFloat(0, 60000.0) - 30000.0);

			males[i].f.x = (int)(RandomFloat(0, 60000.0) - 30000.0);
			males[i].f.y = (int)(RandomFloat(0, 60000.0) - 30000.0);

			females[i].m.x = (int)(RandomFloat(0, 60000.0) - 30000.0);
			females[i].m.y = (int)(RandomFloat(0, 60000.0) - 30000.0);

			females[i].f.x = (int)(RandomFloat(0, 60000.0) - 30000.0);
			females[i].f.y = (int)(RandomFloat(0, 60000.0) - 30000.0);
			//population[i].female.x = RandomFloat(0, 6.0) + surface_min_x;
			//population[i].female.y = RandomFloat(0, 6.0) + surface_min_y;
		}
	}




	void ComputeFitness()
	{
		float max_fitness_males = 0.0f;
		float max_fitness_females = 0.0f;
		sumFitness_males = 0.0f;
		sumFitness_females = 0.0f;
		for (int i = 0; i < males.size(); i++)
		{
			males[i].fitness = GeneHeight_Males(males[i].m) + GeneHeight_Females(males[i].f);
			females[i].fitness = GeneHeight_Males(females[i].m) + GeneHeight_Females(females[i].f);

			if (!GeneInBounds(males[i].m))
			{
				males[i].fitness = 0.0f;
			}
			if (!GeneInBounds(females[i].f))
			{
				females[i].fitness = 0.0f;
			}

			if (males[i].fitness > max_fitness_males)
			{
				max_fitness_males = males[i].fitness;
				best_index_males = i;
			}

			if (females[i].fitness > max_fitness_females)
			{
				max_fitness_females = females[i].fitness;
				best_index_females = i;
			}



			sumFitness_males += males[i].fitness;
			sumFitness_females += females[i].fitness;
		}
	}

	int ChooseParent(int parent_to_skip, vector<individual>& population, float sumFitness, int best_index)
	{
		int randSelector = (int)RandomFloat(0, sumFitness) * 100;

		int memberID = 0;
		int partialSum = 0;

		for (int i = 0; i < population.size(); i++) {

			int n = (int)(population[i].fitness * 100);
			if (partialSum < randSelector && partialSum + n >= randSelector)
			{
				if (i == parent_to_skip)
				{
					if (i + 1 == population.size()) return best_index; // more breeding with population best
					else return i + 1;
				}
				else
				{
					return i;
				}
			}
			partialSum += n;
		}

		return best_index;
	}


	void Breed()
	{
		CopyPopulation(males, mating_pool_males);
		CopyPopulation(females, mating_pool_females);
		bool same_parent = true;
		for (int i = 0; i < males.size(); i++)
		{
			if (!same_parent)
			{
				if (i != best_index_males)
				{
					int m_a = ChooseParent(i, mating_pool_males, sumFitness_males, best_index_males);
					int m_b = ChooseParent(-1, mating_pool_females, sumFitness_females, best_index_females);
					GeneCrossOverSinglePoint(males[i].m, mating_pool_males[m_a].m, mating_pool_females[m_b].m);
					GeneCrossOverSinglePoint(males[i].f, mating_pool_males[m_a].f, mating_pool_females[m_b].f);

				}

				if (i != best_index_females)
				{
					int f_a = ChooseParent(-1, mating_pool_males, sumFitness_males, best_index_males);
					int f_b = ChooseParent(i, mating_pool_females, sumFitness_females, best_index_females);
					GeneCrossOverSinglePoint(females[i].m, mating_pool_males[f_a].m, mating_pool_females[f_b].m);
					GeneCrossOverSinglePoint(females[i].f, mating_pool_males[f_a].f, mating_pool_females[f_b].f);
				}
			}
			else
			{
				// the male and female parents make 1 son and 1 daughter

				int m_a = ChooseParent(i, mating_pool_males, sumFitness_males, best_index_males);
				int m_b = ChooseParent(i, mating_pool_females, sumFitness_females, best_index_females);

				if (i != best_index_males)
				{

					GeneCrossOverSinglePoint(males[i].m, mating_pool_males[m_a].m, mating_pool_females[m_b].m);
					GeneCrossOverSinglePoint(males[i].f, mating_pool_males[m_a].f, mating_pool_females[m_b].f);

				}

				if (i != best_index_females)
				{

					GeneCrossOverSinglePoint(females[i].m, mating_pool_males[m_a].m, mating_pool_females[m_b].m);
					GeneCrossOverSinglePoint(females[i].f, mating_pool_males[m_a].f, mating_pool_females[m_b].f);
				}
			}

			MutatePopulation();
		}
	}

	void  UpdateSimulation()
	{
		ComputeFitness();
		Breed();

		//cout << "Best Index Males: " << best_index_males << ", Best Index Females: " << best_index_females << endl;
		//	cout << "H0: " << males[best_index_males].fitness << ", H1: " << females[best_index_females].fitness << endl;

		//	cout << "=============================" << endl;
	}

	void PrintBest()
	{
		cout << "=============================" << endl;
		cout << " SIMULATION COMPLETE" << endl;
		cout << "=============================" << endl;
		float H = PEAKS_SURFACE::computePeaks((float)males[best_index_males].m.x / 10000, (float)males[best_index_males].m.y / 10000);

		cout << "Best Index: " << best_index_males << ", fitness: " << males[best_index_males].fitness << endl;
		cout << "x: " << (float)males[best_index_males].m.x / 10000 << ", y: " << (float)males[best_index_males].m.y / 10000 << endl;

		DNA_OUTPUT::print_dna(males[best_index_males].m.x);
		cout << ", ";
		DNA_OUTPUT::print_dna(males[best_index_males].m.y);
		cout << endl;

		cout << "NOW RUNNING THE BEST MALES GENES THROUGH THE FEMALE FUNCTION" << endl;
		H = GeneHeight_Females(males[best_index_males].f);
		cout << "H: " << H << ", x: " << (float)males[best_index_males].f.x / 10000 << ", y: " << (float)males[best_index_males].f.y / 10000 << endl;
		cout << "=============================" << endl;





		cout << "Best Index: " << best_index_females << ", fitness: " << females[best_index_females].fitness << endl;
		cout << "x: " << (float)females[best_index_females].f.x / 10000 << ", y: " << (float)females[best_index_females].f.y / 10000 << endl;

		DNA_OUTPUT::print_dna(females[best_index_females].f.x);
		cout << ", ";
		DNA_OUTPUT::print_dna(females[best_index_females].f.y);
		cout << endl;

		cout << "NOW RUNNING THE BEST FEMALES GENES THROUGH THE MALE FUNCTION" << endl;
		H = GeneHeight_Males(females[best_index_females].m);
		cout << "H: " << H << ", x: " << (float)females[best_index_females].m.x / 10000 << ", y: " << (float)females[best_index_females].m.y / 10000 << endl;

		cout << "=============================" << endl;
		cout << " SIMULATION COMPLETE" << endl;
		cout << "=============================" << endl;
	}

	float BestIndexGradient_X()
	{
		return PEAKS_SURFACE::peaksDerivX((float)males[best_index_males].m.x / 10000, (float)males[best_index_males].m.y / 10000);
	}

	float BestIndexGradient_Y()
	{
		return PEAKS_SURFACE::peaksDerivY((float)males[best_index_males].m.x / 10000, (float)males[best_index_males].m.y / 10000);
	}



	bool SIM_STOP()
	{
		//if (population[best_index].fitness < 8.5) return false;

		return false;
	}

}
int _tmain(int argc, _TCHAR* argv[])
{
	srand(5436);

	SURFACE_CLIMBING_GA_2::CreatePopulation(100);
	SURFACE_CLIMBING_GA_2::RamdomizeInitialPopulation();
	for (int epoch = 0; epoch < 1000; epoch++)
	{
		//cout << "Epoch: " << epoch << endl;
		SURFACE_CLIMBING_GA_2::UpdateSimulation();

		if (SURFACE_CLIMBING_GA_2::SIM_STOP())
		{
			SURFACE_CLIMBING_GA_2::PrintBest();
			break;
		}
	}
	SURFACE_CLIMBING_GA_2::PrintBest();
	//SURFACE_CLIMBING_GA::MoveBestToTop();

	// Grid search to prove the GA was successful in 700 epochs
	/*float H = 0.0f;
	float bX = 0.0f, bY = 0.0f;
	for (int i = 0; i < 6000; i++)
	{
	for (int j = 0; j < 6000; j++)
	{
	float X = ((float)((i)  - 3000)) / 1000;
	float Y = ((float)((j) - 3000)) / 1000;
	float f = PEAKS_SURFACE::computePeaks(X, Y);
	if (f > H)
	{
	H = f;
	bX = X;
	bY = Y;
	}
	}
	}

	cout << "bX: " << bX << ", bY: " << bY << ", H: " << H << endl;*/
	return 0;
}

