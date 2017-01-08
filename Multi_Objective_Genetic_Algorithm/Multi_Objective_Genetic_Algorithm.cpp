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
		gene shared;

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

	float GeneHeight_Shared(gene& g)
	{
		float x = ((float)(g.x)) / 10000;
		float y = ((float)(g.y)) / 10000;

		return PEAKS_SURFACE::computePeaks(-x, -y);
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
	float same_parent_probability = 0.5f;


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
				GeneMutate(males[i].shared);
			}
			if (i != best_index_females)
			{
				GeneMutate(females[i].m);
				GeneMutate(females[i].f);
				GeneMutate(females[i].shared);
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

			males[i].shared.x = (int)(RandomFloat(0, 60000.0) - 30000.0);
			males[i].shared.y = (int)(RandomFloat(0, 60000.0) - 30000.0);

			females[i].m.x = (int)(RandomFloat(0, 60000.0) - 30000.0);
			females[i].m.y = (int)(RandomFloat(0, 60000.0) - 30000.0);

			females[i].f.x = (int)(RandomFloat(0, 60000.0) - 30000.0);
			females[i].f.y = (int)(RandomFloat(0, 60000.0) - 30000.0);

			females[i].shared.x = (int)(RandomFloat(0, 60000.0) - 30000.0);
			females[i].shared.y = (int)(RandomFloat(0, 60000.0) - 30000.0);
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
			males[i].fitness = GeneHeight_Males(males[i].m) + GeneHeight_Shared(males[i].shared);
			females[i].fitness = GeneHeight_Females(females[i].f) + GeneHeight_Shared(females[i].shared);

			if (!GeneInBounds(males[i].m) || !GeneInBounds(males[i].shared))
			{
				males[i].fitness = 0.0f;
			}
			if (!GeneInBounds(females[i].f) || !GeneInBounds(females[i].shared))
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
		//bool same_parent = false;
		
		for (int i = 0; i < males.size(); i++)
		{
			float same_parent = 0.6;// RandomFloat(0.0, 1.0);
			if (same_parent > same_parent_probability)
			{
				if (i != best_index_males)
				{
					int m_a = ChooseParent(i, mating_pool_males, sumFitness_males, best_index_males);
					int m_b = ChooseParent(-1, mating_pool_females, sumFitness_females, best_index_females);
					GeneCrossOverSinglePoint(males[i].m, mating_pool_males[m_a].m, mating_pool_females[m_b].m);
					GeneCrossOverSinglePoint(males[i].f, mating_pool_males[m_a].f, mating_pool_females[m_b].f);
					GeneCrossOverSinglePoint(males[i].shared, mating_pool_males[m_a].shared, mating_pool_females[m_b].shared);
				}

				if (i != best_index_females)
				{
					int f_a = ChooseParent(-1, mating_pool_males, sumFitness_males, best_index_males);
					int f_b = ChooseParent(i, mating_pool_females, sumFitness_females, best_index_females);
					GeneCrossOverSinglePoint(females[i].m, mating_pool_males[f_a].m, mating_pool_females[f_b].m);
					GeneCrossOverSinglePoint(females[i].f, mating_pool_males[f_a].f, mating_pool_females[f_b].f);
					GeneCrossOverSinglePoint(females[i].shared, mating_pool_males[f_a].shared, mating_pool_females[f_b].shared);
				}
			}
			else
			{
				// the male and female parents make 1 son and 1 daughter

				int m_a = ChooseParent(i, mating_pool_males, sumFitness_males, best_index_males);
				int m_b = ChooseParent(i, mating_pool_females, sumFitness_females, best_index_females);

				//if (i != best_index_males)
				{

					GeneCrossOverSinglePoint(males[i].m, mating_pool_males[m_a].m, mating_pool_females[m_b].m);
					GeneCrossOverSinglePoint(males[i].f, mating_pool_males[m_a].f, mating_pool_females[m_b].f);
					GeneCrossOverSinglePoint(males[i].shared, mating_pool_males[m_a].shared, mating_pool_females[m_b].shared);
				}

				//if (i != best_index_females)
				{

					GeneCrossOverSinglePoint(females[i].m, mating_pool_males[m_a].m, mating_pool_females[m_b].m);
					GeneCrossOverSinglePoint(females[i].f, mating_pool_males[m_a].f, mating_pool_females[m_b].f);
					GeneCrossOverSinglePoint(females[i].shared, mating_pool_males[m_a].shared, mating_pool_females[m_b].shared);
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

		cout << "Best Index: " << best_index_males << ", H: " << H << endl;
		cout << "x: " << (float)males[best_index_males].m.x / 10000 << ", y: " << (float)males[best_index_males].m.y / 10000 << endl;

		DNA_OUTPUT::print_dna(males[best_index_males].m.x);
		cout << ", ";
		DNA_OUTPUT::print_dna(males[best_index_males].m.y);
		cout << endl;

		cout << "NOW RUNNING THE BEST MALES GENES THROUGH THE FEMALE FUNCTION" << endl;
		H = GeneHeight_Females(males[best_index_males].f);
		cout << "H: " << H << ", x: " << (float)males[best_index_males].f.x / 10000 << ", y: " << (float)males[best_index_males].f.y / 10000 << endl;
		
		cout << "NOW RUNNING THE BEST MALES GENES THROUGH THE SHARED" << endl;
		H = GeneHeight_Shared(males[best_index_males].shared);
		cout << "H: " << H << ", x: " << (float)males[best_index_males].shared.x / 10000 << ", y: " << (float)males[best_index_males].shared.y / 10000 << endl;

		
		
		cout << "=============================" << endl;



		H = GeneHeight_Females(females[best_index_females].f);

		cout << "Best Index: " << best_index_females << ", fitness: " << H << endl;
		cout << "x: " << (float)females[best_index_females].f.x / 10000 << ", y: " << (float)females[best_index_females].f.y / 10000 << endl;

		DNA_OUTPUT::print_dna(females[best_index_females].f.x);
		cout << ", ";
		DNA_OUTPUT::print_dna(females[best_index_females].f.y);
		cout << endl;

		cout << "NOW RUNNING THE BEST FEMALES GENES THROUGH THE MALE FUNCTION" << endl;
		H = GeneHeight_Males(females[best_index_females].m);
		cout << "H: " << H << ", x: " << (float)females[best_index_females].m.x / 10000 << ", y: " << (float)females[best_index_females].m.y / 10000 << endl;

		cout << "NOW RUNNING THE BEST FEMALES GENES THROUGH THE SHARED FUNCTION" << endl;
		H = GeneHeight_Shared(females[best_index_females].shared);
		cout << "H: " << H << ", x: " << (float)females[best_index_females].shared.x / 10000 << ", y: " << (float)females[best_index_females].shared.y / 10000 << endl;

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


namespace GA
{
	using std::vector;
	using namespace BIT_OPS;



	class gene
	{
	public:
		gene(){ }
		
		// accessors
		int GetNumBits() { return 32; }
		int GetGeneAsInt() { return the_int; }
		float GetGeneAsFloat() { return *reinterpret_cast<float*>(&the_int); }
		char GetGeneAsChar(){ return *reinterpret_cast<char*>(&the_int); }
		
		// modifiers
		void SetGene(int X) { the_int = X; }
		void SetGene(float A) { 
			the_int = *reinterpret_cast< int*>(&A);
		}
		void SetGene(char c){
			the_int = *reinterpret_cast< int*>(&c);
		}

		// members
		int the_int;
		int type = 0;
	};



	class GeneOp
	{
	public:


		
		void MutateBit(gene& g, int mut_point)
		{
			int X = 0;
			int num_bits = 32;
			X = g.GetGeneAsInt();

			if (mut_point > num_bits) return;
			
			if (CheckBit(X, mut_point))
			{
				ResetBit(X, mut_point);
			}
			else
			{
				SetBit(X, mut_point);
			}

			g.SetGene(X);
		}

		void CrossOver(gene &g, gene &A, gene &B, int cross_point)
		{
			int num_bits = g.GetNumBits();
			if (cross_point > num_bits) return;

			int X = g.GetGeneAsInt();
			int a = A.GetGeneAsInt();
			int b = B.GetGeneAsInt();

			for (int i = 0; i < cross_point; i++)
			{
				if (CheckBit(a, i))
				{
					SetBit(X, i);
				}
			}
			for (int i = cross_point; i < num_bits; i++)
			{
				if (CheckBit(b, i))
				{
					SetBit(X, i);
				}
			}

			g.SetGene(X);

		}
	};

	


	class Chromosome
	{
	public:
		bool IsCompatible(Chromosome &A, Chromosome &B){ return true; }
		void CrossOverAssignment(Chromosome &A, Chromosome &B, int cross_point)
		{
			for (int i = 0; i < cross_point; i++)
			{
				memcpy((void*)&dna[i], (void*)&A.dna[i], sizeof(gene));
			}
			for (int i = cross_point; i < dna.size(); i++)
			{
				memcpy((void*)&dna[i], (void*)&B.dna[i], sizeof(gene));
			}
		}
		void CrossOverInsert(Chromosome &A, Chromosome &B, int cross_point)
		{
			for (int i = 0; i < cross_point; i++)
			{
				dna.push_back(A.dna[i]);
			}
			for (int i = cross_point; i < dna.size(); i++)
			{
				dna.push_back(B.dna[i]);
			}
		}
		void CrossOver(Chromosome &A, Chromosome &B, int cross_point)
		{
			if (cross_point > A.dna.size() || cross_point > B.dna.size())
				return;

			if ((A.dna.size() == B.dna.size()) ||
				((this->dna.size() < A.dna.size()) && 
				 (this->dna.size() < B.dna.size()) && 
				 (cross_point < this->dna.size())))
			{
				CrossOverAssignment(A, B, cross_point);
			}
			else
			{
				dna.clear();
				CrossOverInsert(A, B, cross_point);
			}

		}

		vector< gene > dna;
	};


	class Individual
	{
	public:
		Chromosome dna;
	};
}


namespace BIT_BUFFER
{
	using namespace DNA_OUTPUT;

	class BlockIndexBitField
	{
	public:

		BlockIndexBitField(int numBits)
		{
			dataSize = 32;
			if (numBits == 0)
				numBits = dataSize;

			m_numBits = numBits;
			m_sizeInInts = (numBits / dataSize);// +1; // need to add one int
			if (numBits % dataSize != 0) m_sizeInInts++;

			m_data = new int[(m_sizeInInts)];

			//memset( (void*)m_data, 0, sizeof(m_data) );
			for (int i = 0; i < m_sizeInInts; i++)
			{
				m_data[i] = 0;
			}
		}

		~BlockIndexBitField()
		{
			if (m_data)
				delete[] m_data;
		}

		void SetBit(int idx)
		{
			if (idx > m_numBits) return;

			if (idx < 0) return;

			int blockIndex = (int)((float)(idx) / dataSize);

			int remainder = idx - blockIndex*dataSize;

			m_data[blockIndex] |= (1 << remainder);
		}

		bool CheckBit(int idx)
		{
			if (idx > m_numBits) return false;

			if (idx < 0) return false;

			int blockIndex = (int)((float)(idx) / dataSize);

			int remainder = idx - blockIndex*dataSize;

			return m_data[blockIndex] & (1 << remainder);
		}

		void ResetBit(int idx)
		{
			if (idx > m_numBits) return;

			if (idx < 0) return;

			int blockIndex = (int)((float)(idx) / dataSize);

			int remainder = idx - blockIndex*dataSize;

			m_data[blockIndex] &= ~(1 << remainder);
		}

		void Clear()
		{
			for (int i = 0; i < m_sizeInInts; i++)
			{
				m_data[i] = 0;
			}
		}

		bool Binary_AND(BlockIndexBitField &B, BlockIndexBitField &C)
		{
			if (B.m_numBits != C.m_numBits)
			{
				return false;
			}

			if (this->m_numBits != B.m_numBits)
			{
				return false;
			}

			for (int i = 0; i < B.m_sizeInInts; i++)
			{
				this->m_data[i] = B.m_data[i] & C.m_data[i];
			}

			return true;
		}

		bool Binary_OR(BlockIndexBitField &B, BlockIndexBitField &C)
		{
			if (B.m_numBits != C.m_numBits)
			{
				return false;
			}

			if (this->m_numBits != B.m_numBits)
			{
				return false;
			}

			for (int i = 0; i < B.m_sizeInInts; i++)
			{
				this->m_data[i] = B.m_data[i] | C.m_data[i];
			}

			return true;
		}

		bool Binary_NAND(BlockIndexBitField &B, BlockIndexBitField &C)
		{
			if (B.m_numBits != C.m_numBits)
			{
				return false;
			}

			if (this->m_numBits != B.m_numBits)
			{
				return false;
			}

			for (int i = 0; i < B.m_sizeInInts; i++)
			{
				this->m_data[i] = B.m_data[i] & ~C.m_data[i];
			}

			return true;
		}

		bool Binary_XOR(BlockIndexBitField &B, BlockIndexBitField &C)
		{
			if (B.m_numBits != C.m_numBits)
			{
				return false;
			}

			if (this->m_numBits != B.m_numBits)
			{
				return false;
			}

			for (int i = 0; i < B.m_sizeInInts; i++)
			{
				this->m_data[i] = B.m_data[i] ^ C.m_data[i];
			}

			return true;
		}

		void print()
		{
			for (int i = 0; i < m_sizeInInts; i++)
			{
				DNA_OUTPUT::print_dna(m_data[i]);
			}
		}

		int m_sizeInInts;
		int m_numBits;
		int dataSize;
		int *m_data;
	};

}

#include <string>
#include <unordered_map>

using std::string;
using std::unordered_map;

namespace MICRO_BIOLOGY
{
	using namespace BIT_OPS;
	using namespace DNA_OUTPUT;
	using namespace BIT_BUFFER;

	enum NUCLEOTIDES { U, C, A, G };

	unordered_map< char, NUCLEOTIDES > m;

	string AMINO_ACIDS[4][4][4];

	int NUM_SYMBOLS_IN_AMINO = 3;
	int NUM_BITS_FOR_NUCLEOTIDE = 2;

	BlockIndexBitField bitfield(32 * NUM_SYMBOLS_IN_AMINO * NUM_BITS_FOR_NUCLEOTIDE);
	
	int bit_sequence[6] = { 0, 0, 0, 0, 0, 0 };

	string s = "";

	void Initialize()
	{
		
		for (int i = 0; i < 96; i++)
		{
			int R = RandomInt(0, 4);
			if (R == 0)
			{
				s += "G";
			}
			if (R == 1)
			{
				s += "A";
			}
			if (R == 2)
			{
				s += "T";
			}
			if (R == 3)
			{
				s += "C";
			}
		}

		m['T'] = U;
		m['C'] = C;
		m['A'] = A;
		m['G'] = G;



		AMINO_ACIDS[U][U][U] = "Phenylalanine";
		AMINO_ACIDS[U][U][C] = "Phenylalanine";
		AMINO_ACIDS[U][U][A] = "Leucine";
		AMINO_ACIDS[U][U][G] = "Leucine";

		AMINO_ACIDS[U][C][U] = "Serine";
		AMINO_ACIDS[U][C][C] = "Serine";
		AMINO_ACIDS[U][C][A] = "Serine";
		AMINO_ACIDS[U][C][G] = "Serine";

		AMINO_ACIDS[U][A][U] = "Phenylalanine";
		AMINO_ACIDS[U][A][C] = "Phenylalanine";
		AMINO_ACIDS[U][A][A] = "stop";
		AMINO_ACIDS[U][A][G] = "stop";

		AMINO_ACIDS[U][G][U] = "Cysteine";
		AMINO_ACIDS[U][G][C] = "Cysteine";
		AMINO_ACIDS[U][G][A] = "stop";
		AMINO_ACIDS[U][G][G] = "Tryptophan";


		AMINO_ACIDS[C][U][U] = "Leucine";
		AMINO_ACIDS[C][U][C] = "Leucine";
		AMINO_ACIDS[C][U][A] = "Leucine";
		AMINO_ACIDS[C][U][G] = "Leucine";

		AMINO_ACIDS[C][C][U] = "Proline";
		AMINO_ACIDS[C][C][C] = "Proline";
		AMINO_ACIDS[C][C][A] = "Proline";
		AMINO_ACIDS[C][C][G] = "Proline";

		AMINO_ACIDS[C][A][U] = "Histidine";
		AMINO_ACIDS[C][A][C] = "Histidine";
		AMINO_ACIDS[C][A][A] = "Glutamine";
		AMINO_ACIDS[C][A][G] = "Glutamine";

		AMINO_ACIDS[C][G][U] = "Arginine";
		AMINO_ACIDS[C][G][C] = "Arginine";
		AMINO_ACIDS[C][G][A] = "Arginine";
		AMINO_ACIDS[C][G][G] = "Arginine";


		AMINO_ACIDS[A][U][U] = "Isolecine";
		AMINO_ACIDS[A][U][C] = "Isolecine";
		AMINO_ACIDS[A][U][A] = "Isolecine";
		AMINO_ACIDS[A][U][G] = "Methionine";

		AMINO_ACIDS[A][C][U] = "Threonine";
		AMINO_ACIDS[A][C][C] = "Threonine";
		AMINO_ACIDS[A][C][A] = "Threonine";
		AMINO_ACIDS[A][C][G] = "Threonine";

		AMINO_ACIDS[A][A][U] = "Asparagine";
		AMINO_ACIDS[A][A][C] = "Asparagine";
		AMINO_ACIDS[A][A][A] = "Lysine";
		AMINO_ACIDS[A][A][G] = "Lysine";

		AMINO_ACIDS[A][G][U] = "Serine";
		AMINO_ACIDS[A][G][C] = "Serine";
		AMINO_ACIDS[A][G][A] = "Arginine";
		AMINO_ACIDS[A][G][G] = "Arginine";

		AMINO_ACIDS[G][U][U] = "Valine";
		AMINO_ACIDS[G][U][C] = "Valine";
		AMINO_ACIDS[G][U][A] = "Valine";
		AMINO_ACIDS[G][U][G] = "Valine";

		AMINO_ACIDS[G][C][U] = "Alanine";
		AMINO_ACIDS[G][C][C] = "Alanine";
		AMINO_ACIDS[G][C][A] = "Alanine";
		AMINO_ACIDS[G][C][G] = "Alanine";

		AMINO_ACIDS[G][A][U] = "Aspartate";
		AMINO_ACIDS[G][A][C] = "Aspartate";
		AMINO_ACIDS[G][A][A] = "Glutamate";
		AMINO_ACIDS[G][A][G] = "Glutamate";

		AMINO_ACIDS[G][G][U] = "Glycine";
		AMINO_ACIDS[G][G][C] = "Glycine";
		AMINO_ACIDS[G][G][A] = "Glycine";
		AMINO_ACIDS[G][G][G] = "Glycine";
	}

	void TestDNAString(string s)
	{

		cout << s << endl;
		
		int idx = 0;

		for (int i = 0; i < s.length(); i++)
		{
			idx = i / 16;
			int b = 2*(i - 16 * idx);

			if (s[i] == 'G')
			{
				cout << "00";
				
			}
			if (s[i] == 'A')
			{
				cout << "01"; 
				SetBit(bit_sequence[idx], b + 1);
			}
			if (s[i] == 'T')
			{
				cout << "10";
				SetBit(bit_sequence[idx], b);
			}
			if (s[i] == 'C')
			{
				cout << "11";
				SetBit(bit_sequence[idx], b);
				SetBit(bit_sequence[idx], b + 1);
			}
		}


		cout << endl;
		for (int i = 0; i < 6; DNA_OUTPUT::print_dna(bit_sequence[i++]));
		cout << endl;

		std::vector< std::string > aminos;
		std::string curr_amino = "";
		for (int i = 0; i < 96; i++)
		{
				idx = i / 16;
				int b =  2*(i - 16 * idx);

				if (CheckBit(bit_sequence[idx], b))
				{
					if (CheckBit(bit_sequence[idx], b + 1))
					{
						curr_amino += "C";
					}
					else
					{
						curr_amino += "T";
					}
				}
				else
				{
					if (CheckBit(bit_sequence[idx], b + 1))
					{
						curr_amino += "A";
					}
					else
					{
						curr_amino += "G";
					}
				}

				if (curr_amino.length() == 3)
				{

					aminos.push_back(curr_amino);
					curr_amino = "";
				}
			
		}

		cout << aminos.size() << endl;


		cout << endl;

		std::string am = "";
		for (int i = 0; i < aminos.size(); i++)
		{
			cout << "code: " << aminos[i];

			int idx1 = (int)m[ aminos[i][0] ];
			int idx2 = (int)m[ aminos[i][1] ];
			int idx3 = (int)m[ aminos[i][2] ];
			cout << ",   amino acid: " << AMINO_ACIDS[idx1][idx2][idx3] << endl;

			am += aminos[i];

		}

		cout << endl;
		if (am != s) cout << "am is not equal to s " << endl;
	}


	void TestDNAstring_BlockIndexBitfield(string s)
	{
	
		cout << s << endl;

		int idx = 0;

		for (int i = 0; i < s.length(); i++)
		{
			int b = 2 * i;
			//cout << idx << b << " ";

			if (s[i] == 'G')
			{
				cout << "00";

			}
			if (s[i] == 'A')
			{
				cout << "01";

				bitfield.SetBit(b + 1);
			}
			if (s[i] == 'T')
			{
				cout << "10";
				bitfield.SetBit(b);
			}
			if (s[i] == 'C')
			{
				cout << "11";
				bitfield.SetBit(b);
				bitfield.SetBit(b + 1);
			}
		}

		
		cout << endl;
		bitfield.print();
		cout << endl;

		std::vector< std::string > aminos;
		std::string curr_amino = "";
		for (int i = 0; i < 96; i++)
		{
			int b = 2 * i;

			if (bitfield.CheckBit( b))
			{
				if (bitfield.CheckBit(b+1))
				{
					curr_amino += "C";
				}
				else
				{
					curr_amino += "T";
				}
			}
			else
			{
				if (bitfield.CheckBit(b + 1))
				{
					curr_amino += "A";
				}
				else
				{
					curr_amino += "G";
				}
			}

			if (curr_amino.length() == 3)
			{

				aminos.push_back(curr_amino);
				curr_amino = "";
			}
		}

		cout << aminos.size() << endl;


		cout << endl;

		std::string am = "";

		for (int i = 0; i < aminos.size(); i++)
		{
			cout << "code: " << aminos[i];

			int idx1 = (int)m[aminos[i][0]];
			int idx2 = (int)m[aminos[i][1]];
			int idx3 = (int)m[aminos[i][2]];
			cout << ",   amino acid: " << AMINO_ACIDS[idx1][idx2][idx3] << endl;
			am += aminos[i];

		}

		cout << endl;
		if (am != s) cout << "am is not equal to s " << endl;
	}
}

#include <list>

using std::list;


int _tmain(int argc, _TCHAR* argv[])
{

	MICRO_BIOLOGY::Initialize();
	MICRO_BIOLOGY::TestDNAstring_BlockIndexBitfield(MICRO_BIOLOGY::s);
	MICRO_BIOLOGY::TestDNAString(MICRO_BIOLOGY::s);

/*	There might be a reason to use a list instead of a vector 
    if I change to large genomes with growable genes 
	http://baptiste-wicht.com/posts/2012/12/cpp-benchmark-vector-list-deque.html

list< int > a_list;8

	a_list.push_back(1);
	a_list.push_back(2);
	a_list.push_back(8);

	for (int i : a_list) { std::cout << i << endl; }

	auto i1 = a_list.begin();
	*/

	// QUICK TEST OF NEW GENETIC OBJECTS
	GA::Chromosome c;
	GA::gene g1,g2,g3;

	g1.type=0;
	g2.type = 1;
	g3.type = 2;

	g1.SetGene(25);
	g2.SetGene(1.8f);
	g3.SetGene((int)'a');

	c.dna.push_back(g1);
	c.dna.push_back(g2);
	c.dna.push_back(g3);
	c.dna.push_back(g1);
	c.dna.push_back(g2);
	c.dna.push_back(g2);

	/*for (int i = 0; i < c.dna.size()/2; i++)
	{
		switch (c.dna[i].type)
		{
		case 0: ((GA::int_gene*)&c.dna[i])->SetGene(RandomInt(0, 10)); break;
		case 1:((GA::float_gene*)&c.dna[i])->SetGene(RandomInt(0, 10)); break;
		case 2:((GA::char_gene*)&c.dna[i])->SetGene(RandomInt(0, 10)); break;
		}
	}*/
	for (int i = 0; i < c.dna.size(); i++)
	{
		switch (c.dna[i].type)
		{
		case 0:std::cout << c.dna[i].GetGeneAsInt() << endl; break;
		case 1:{
			int X = c.dna[i].GetGeneAsInt();
			std::cout << *reinterpret_cast<float*>(&X) << endl; break;
		}
		case 2:std::cout << (char)c.dna[i].GetGeneAsInt() << endl; break;
		}
		
	}
	


	// TEST OF HILL CLIMBING GA.
	srand(5436);

	SURFACE_CLIMBING_GA_2::CreatePopulation(100);
	SURFACE_CLIMBING_GA_2::RamdomizeInitialPopulation();
	for (int epoch = 0; epoch < 7000; epoch++)
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

