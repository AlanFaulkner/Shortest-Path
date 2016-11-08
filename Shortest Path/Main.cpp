#include <iostream>
#include <iomanip>
#include <random>
#include <algorithm>
#include <vector>
#include <ctime>

struct Dynamic_Sort {
	Dynamic_Sort(int paramA) { this->paramA = paramA; }
	bool operator () (std::vector<double> i, std::vector<double> j) { return i[paramA] < j[paramA]; }

	int paramA;
};

bool Is_Present(std::vector<std::vector<double>> Results, std::vector<double> Permintation, int Size) {
	
	for (int i = 0; i < Results.size(); i++) {
		std::vector<double> test;
		for (int j = 0; j < Size; j++) {
			test.push_back(Results[i][j]);
		}
		if (std::equal(Permintation.begin(), Permintation.end(), test.begin())) { return true; }
		std::reverse(Permintation.begin(),Permintation.end());
		if (std::equal(Permintation.begin(), Permintation.end(), test.begin())) { return true; }
	}
	return false;
}

bool Double_Sort(std::vector<double>i, std::vector<double>j) {
	if (i[0] < j[0]) { return true; }
	else if (i[0]>j[0]){ return false; }
	else if (i[3] > j[3]) { return true; }
	return false;
}

//below is what is currently in use

std::vector<std::vector<double>> Cities = { { 0,1,4 },{ 1,3,3 },{ 2,1,7 },{ 3,8,9 },{ 4,5,6 },{ 5,9,0 },{ 6,1,1 },{ 7,0,0 },{ 8,6,2 },{ 9,2,2 },
{ 10,1,7,2 },{ 11,3,-3,-3 },{ 12,1,11,2 },{ 13,-8,9,-7 },{ 14,5,-6,6 },{ 15,9,0,-3 },{ 16,-1,-1,0 },{ 17,9,0,-1 },{ 18,6,2,8 },{ 19,-2,2,2 } };

double Get_Path_Length(std::vector<std::vector<double>> Input_Coordinates, int i, int j)
{
	return sqrt(((Input_Coordinates[i][1] - Input_Coordinates[j][1])*(Input_Coordinates[i][1] - Input_Coordinates[j][1]))
		+ ((Input_Coordinates[i][2] - Input_Coordinates[j][2])*(Input_Coordinates[i][2] - Input_Coordinates[j][2])));
}

void Brute_Force_Methood(std::vector<std::vector<double>> Cities, int Number_Of_Cities) {

	/*
	##########################
	## Brute Force Approach ##
	##########################
	*/

		/*
		******************
		** Build inputs **
		******************
		*/
			std::vector<double> Journey;
			for (int i = 0; i < Number_Of_Cities; i++) { Journey.push_back(i); }

		/*
		****************
		** Main loop  **
		****************
		*/

			std::vector<std::vector<double>> Results;

			do {

				double Length = 0;
				//need a tempory vector to store current permutation
				std::vector<double>Current_Journey;
				for (int i = 0; i < Number_Of_Cities; i++) { Current_Journey.push_back(Journey[i]); }
				for (int i = 1; i < Number_Of_Cities; i++) {
					Length += Get_Path_Length(Cities, Current_Journey[i], Current_Journey[i - 1]);
				}
				Current_Journey.push_back(Length);
				Results.push_back(Current_Journey);
			} while (std::next_permutation(Journey.begin(), Journey.end()));

			//Order Results based on path length - shortest first
			std::sort(Results.begin(), Results.end(), Dynamic_Sort(Number_Of_Cities));

		/*
		*******************
		** Print Results **
		*******************
		*/

			std::cout << "##########################" << std::endl << "## Brute Force Approach ##" << std::endl << "##########################" << std::endl << std::endl;
			std::cout << "The shortest path length is: " << Results[0][Number_Of_Cities] << std::endl;
			std::cout << "The Journey taken was: ";
			for (int i = 0; i < Number_Of_Cities; i++) { std::cout << Results[0][i] << " "; }

			return;
}

std::vector<double> Nearest_Neighbour(std::vector < std::vector<double >> Data, int Number_Of_Cities, int Starting_City) {
	
	/*
	################################
	## Nearest Neighbour Approach ##
	################################
	*/

		/*
		******************
		** Build inputs **
		******************
		*/

			std::vector<double> List_Of_Cities;
			for (int i = 0; i < Number_Of_Cities; i++) { List_Of_Cities.push_back(i); }

			std::vector<std::vector<double>>Results;
			std::vector<double>Results_Row;
			std::vector<std::vector<double>> Neighbourhood;
			std::vector<double>Dist_Neighbour;
			std::vector<double> Journey_Data; // Information about path taken and its length.
		/*
		***************
		** Main Loop **
		***************
		*/

			for (int i = 0; i < Number_Of_Cities - 1; i++) {

				// Move starting city to end of vector and delete it from list
				std::remove(List_Of_Cities.begin(), List_Of_Cities.end(), Starting_City);
				List_Of_Cities.pop_back();

				// Get the distance between starting city and each of its neighbours
				for (int i = 0; i < List_Of_Cities.size(); i++) {
					double sum = Get_Path_Length(Data, List_Of_Cities[i], Starting_City);
					Dist_Neighbour.push_back(Starting_City);
					Dist_Neighbour.push_back(List_Of_Cities[i]);
					Dist_Neighbour.push_back(sum);
					Neighbourhood.push_back(Dist_Neighbour);
					Dist_Neighbour.clear();
				}

				// Order cities based on distance from starting city. closest vs farthest
				std::sort(Neighbourhood.begin(), Neighbourhood.end(), Dynamic_Sort(2));

				//Put the best distance into final results table.
				for (int i = 0; i < 3; i++) { Results_Row.push_back(Neighbourhood[0][i]); }
				Results.push_back(Results_Row);

				Starting_City = Results[i][1];
				Neighbourhood.clear();
				Results_Row.clear();
			}

			// Get final path legnth
			double Path_length = 0;
			for (int i = 0; i < Results.size(); i++) { Path_length += Results[i][2]; }
			
			// Build vector containing info about journey.
			Journey_Data.push_back(Results[0][0]);
			for (int i = 0; i < Results.size(); i++) { Journey_Data.push_back(Results[i][1]); }
			Journey_Data.push_back(Path_length);

	return Journey_Data;
}

void Genetic_Algorithm(std::vector<std::vector<double>> Data, int Gene_Size, int Gene_Pool_Size, int Mutation_Rate) {
	
	/*
	#######################
	## Genetic Algorithm ##
	#######################
	*/

	std::vector<std::vector<double>> Gene_Pool;
	std::vector<double> Gene;
	std::vector<std::vector<double>> New_Generation;
	std::vector<double> Child;
	int Condition = 0;
	int Steps = 0;
	/*
	*********************
	** Build Gene Pool **
	*********************
	*/

	for (int Num_Genes = 0; Num_Genes < Gene_Pool_Size; Num_Genes++) {
		// build indiviual genes (note all genes start with 0)
		for (int i = 0; i < Gene_Size; i++) { Gene.push_back(i); }
		std::shuffle(Gene.begin() + 1, Gene.end(), std::random_device());
		//check that new gene is unique
		while (Is_Present(Gene_Pool, Gene, Gene_Size) && Num_Genes>0) {
			std::shuffle(Gene.begin() + 1, Gene.end(), std::random_device());
		}
		//get path length of gene
		double sum = 0;
		for (int i = 1; i < Gene.size(); i++) {
			sum += Get_Path_Length(Data, Gene[i], Gene[i - 1]);
		}
		Gene.push_back(sum);
		Gene_Pool.push_back(Gene);
		Gene.clear();
	}

	//Sort Genes according to 
	std::sort(Gene_Pool.begin(), Gene_Pool.end(), Dynamic_Sort(Gene_Size));

	/*
	***************
	** Main Loop **
	***************
	*/

	//Elvolve new genes

	do {
		std::random_device rd;
		std::mt19937 generator(rd());
		std::uniform_int_distribution<int> distribution(0, 5);
		std::uniform_int_distribution<int> dis(0, 10);
		std::uniform_int_distribution<int> Gene(1, Gene_Size - 1);

		//eliteism - copy best gene from previous gen to new gen
		New_Generation.push_back(Gene_Pool[0]);

		for (int i = 1; i < Gene_Pool_Size; i++) {
			int Mother = distribution(generator);
			int Farther = dis(generator);

			//make sure parents are different;
			while (Farther == Mother) { Farther = distribution(generator); }

			//crossover - takes first first half of the genes from mother then adds the farthers genes to the end in the order they appear in dad
			for (int i = 0; i < Gene_Size%2; i++) { Child.push_back(Gene_Pool[Mother][i]); }
			for (int i = 0; i < Gene_Size; i++) {
				if (std::find(Child.begin(), Child.end(), Gene_Pool[Farther][i]) == Child.end()) { Child.push_back(Gene_Pool[Farther][i]); }
			}

			//mutation swap any two 'genes'
			if (bool Follow = (rand() % 100) < Mutation_Rate) {
				std::swap(Child[Gene(generator)], Child[Gene(generator)]);
			}

			//calcuate new path length of child
			double sum = 0;
			for (int i = 1; i < Child.size(); i++) { sum += Get_Path_Length(Data, Child[i], Child[i - 1]); }

			Child.push_back(sum);
			New_Generation.push_back(Child);
			Child.clear();
		}

		//update Gene pool
		Gene_Pool = New_Generation;
		std::sort(Gene_Pool.begin(), Gene_Pool.end(), Dynamic_Sort(Gene_Size));
		New_Generation.clear();
		Steps++;

		//termination - gene pool contains no unique genes.
		Condition = 0;
		for (int x = 1; x < Gene_Pool.size(); x++) {
			if (Gene_Pool[x][Gene_Size] != Gene_Pool[0][Gene_Size]) { Condition++; }
		}

	} while (Condition != 0);

	/*
	*******************
	** Print Results **
	*******************
	*/

	std::cout << "########################" << std::endl << "## Genetic Algorithum ##" << std::endl << "########################" << std::endl << std::endl;
	std::cout << "Simulation completed in " << Steps << " steps" << std::endl;
	std::cout << "Final path length = " << Gene_Pool[0][Gene_Size] << std::endl;
	std::cout << "Travel order: ";
	for (int j = 0; j < Gene_Size; j++) { std::cout << " " << Gene_Pool[0][j]; }

	return;
}

void Ant_Colony(std::vector<std::vector<double>>Data, int number_of_cities, int number_of_ants, int Follow_Rate, double Decay_Rate, int Steps, bool Full_Print) {

	/*
	################################
	## Ant colony solution to TSP ##
	################################
	*/

	/*
	***************************
	** Simulation Parameters **
	***************************
	*/
	std::vector<std::vector<double>> Map;
	std::vector<std::vector<double>> Colony;
	std::vector<double> Ant;
	std::vector<double> Map_Data;

	double Current_Shortest_Path = 0;
	double Global_Shortest_Path = 0;
	std::vector<double> Path;

	/*
	******************
	** Build inputs **
	******************
	*/

	// make a map of all cities in problem, the distance between them and set an inital pheromone level to 0
	for (int i = 0; i < number_of_cities; i++) {
		for (int j = 0; j < number_of_cities; j++) {

			//add the two cities
			Map_Data.push_back(i);
			Map_Data.push_back(j);
			//distance between the two cities
			Map_Data.push_back(Get_Path_Length(Data, j, i));
			//inital pharomone level is 0
			Map_Data.push_back(0);
			Map.push_back(Map_Data);
			Map_Data.clear();
		}
	}

	//send each ant in colony on a random trip and save information and total journey length

	for (int i = 0; i < number_of_ants; i++) {
		for (int j = 0; j < number_of_cities; j++) { Ant.push_back(j); }
		std::shuffle(Ant.begin() + 1, Ant.end(), std::random_device());
		double jouney_length = 0;
		for (int x = 1; x < number_of_cities; x++) { jouney_length += Get_Path_Length(Data, Ant[x], Ant[x - 1]); }
		Ant.push_back(jouney_length);
		Colony.push_back(Ant);
		Ant.clear();
	}

	std::sort(Colony.begin(), Colony.end(), Dynamic_Sort(number_of_cities));
	Global_Shortest_Path = Colony[0][number_of_cities];
	for (int i = 0; i < number_of_cities; i++) { Path.push_back(Colony[0][i]); }
	Current_Shortest_Path = Global_Shortest_Path;

	/*
	***************
	** Main Loop **
	***************
	*/

	for (int steps = 0; steps < Steps; steps++) {

		//update pharomone on city map - this is done by dividing food(value 100) by the total journey length of a given ant by the distance
		//between cities

		for (int i = 0; i < number_of_ants; i++) {
			for (double j = 1; j < number_of_cities; j++) {

				int Start_City = Colony[i][j - 1];
				int End_City = (Start_City*(number_of_cities)) + Colony[i][j];
				Map[End_City][3] += (10 / Colony[i][number_of_cities]);

			}
		}

		if (Full_Print == true) {
			//test print
			for (int i = 0; i < Map.size(); i++) {
				for (int j = 0; j < 4; j++) {
					std::cout << std::setprecision(3) << std::setw(5) << Map[i][j] << " ";
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
			for (int i = 0; i < Colony.size(); i++) {
				for (int j = 0; j < number_of_cities + 1; j++) {
					std::cout << std::setprecision(3) << std::setw(5) << Colony[i][j] << " ";
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}

		//Ants now move between cities by choosing which city has the highest value and any other city by random

		//generates a 50:50 probability
		std::random_device rd;
		std::mt19937 generator(rd());
		std::uniform_int_distribution<int> City(0, number_of_cities - 1);
		Colony.clear();
		for (int i = 0; i < number_of_ants; i++) {

			Ant.push_back(0);

			for (int j = 1; j < number_of_cities; j++) {

				if (bool Follow = (rand() % 100) < Follow_Rate) {
					int current_city = Ant[Ant.size() - 1];
					int map_positon = (current_city*(number_of_cities));
					while (std::find(Ant.begin(), Ant.end(), Map[map_positon][1]) != Ant.end()) { map_positon++; }
					Ant.push_back(Map[map_positon][1]);
				}
				else {
					//eliminate any cities already visited and picks next city at random.
					int which_city = City(generator);
					while (std::find(Ant.begin(), Ant.end(), which_city) != Ant.end()) { which_city = City(generator); }
					Ant.push_back(which_city);
				}
			}

			//add path length and add info to ant colony
			double jouney_length = 0;
			for (int x = 1; x < number_of_cities; x++) { jouney_length += Get_Path_Length(Data, Ant[x], Ant[x - 1]); }
			Ant.push_back(jouney_length);
			Colony.push_back(Ant);
			Ant.clear();

			//pharomone exponetial decay rate
			for (int i = 0; i < Map.size(); i++) {
				Map[i][3] = Map[i][3] * Decay_Rate;
				if (Map[i][3] < 0.001) { Map[i][3] = 0; }
			}

		}

		std::sort(Colony.begin(), Colony.end(), Dynamic_Sort(number_of_cities));
		Current_Shortest_Path = Colony[0][number_of_cities];
		if (Current_Shortest_Path < Global_Shortest_Path) {
			Global_Shortest_Path = Current_Shortest_Path;
			for (int i = 0; i < number_of_cities; i++) { Path[i] = Colony[0][i]; }
		}
	}

	/*
	*******************
	** Print Results **
	*******************
	*/

	std::cout << "#########################" << std::endl << "## Ant Colony Solution ##" << std::endl << "#########################" << std::endl << std::endl;
	std::cout << "Simulaton Length: " << Steps << std::endl;
	std::cout << "Colony Size: " << number_of_ants << std::endl;
	std::cout << "Pharomone Decay Rate: " << Decay_Rate << std::endl;
	std::cout << "Follow Rate: " << Follow_Rate << " %"<< std::endl << std::endl;
	std::cout << "The shortest path length is: " << Global_Shortest_Path << std::endl;
	std::cout << "The Journey taken was: ";
	for (int i = 0; i < Path.size(); i++) { std::cout << Path[i] << " "; }

	return;
}
//
//int main() {
//
///*
//######################################
//## Traveling Salesman Problem (TSP) ##
//######################################
//*/
//
//	/*
//	***************************
//	** Simulation Parameters **
//	***************************
//	*/
//
//		int Number_Of_Cities = 5;
//
//	/*
//	##########################
//	## Brute Force Approach ##
//	##########################
//	*/
//
//	Brute_Force_Methood(Cities, Number_Of_Cities);
//	std::cout << std::endl << std::endl;
//
//	/*
//	################################
//	## Nearest Neighbour Approach ##
//	################################
//	*/
//
//	Nearest_Neighbour(Cities, Number_Of_Cities, 3);
//	std::cout << std::endl << std::endl;
//
//	/*
//	#############################
//	## Genetic solution to TSP ##
//	#############################
//	*/
//	
//	Genetic_Algorithm(Cities, Number_Of_Cities, Number_Of_Cities * 4,10);
//	std::cout << std::endl << std::endl;
//
//	/*
//	################################
//	## Ant colony solution to TSP ##
//	################################
//	*/
//	
//	Ant_Colony(Cities, Number_Of_Cities, 40, 70, 0.6, 5000, false);
//
//	return 0;
//}

int main() {
	clock_t start, end;
	for (int i = 2; i < 11; i++) {
		start = clock();
		Brute_Force_Methood(Cities, i);
		end = clock();
		double msecs = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;
		std::cout << std::endl << "Time taken: " << msecs << std::endl << std::endl;
	}
}

//int main() {
//	clock_t start, end;
//	for (int i = 2; i < 11; i++) {
//		start = clock();
//		std::vector<std::vector<double>> Shorest_Path;
//		for (int j = 0; j < i; j++) { // There are i possible number of cities to start at.
//			Shorest_Path.push_back(Nearest_Neighbour(Cities, i, j));
//		}
//		std::sort(Shorest_Path.begin(),Shorest_Path.end(),Dynamic_Sort(i)); // Shorest journey first.
//		end = clock();
//		double msecs = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;
//		std::cout << std::endl << "Shorest Path: " << Shorest_Path[0][i];
//		std::cout << std::endl << "Journey made: ";
//		for (int a = 0; a < i; a++) { std::cout << Shorest_Path[0][a] << " "; }
//		std::cout << std::endl << "Time taken: " << msecs << std::endl << std::endl;
//	}
//}



