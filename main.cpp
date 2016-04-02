#include <iostream>
#include <iomanip>
#include <conio.h>
#include <math.h>
#include <algorithm>
#include <string>
#include <fstream>
#include <map>
using namespace std;

/* CONSTANTS */
const float R = 0.0019872; /* (kcal) */
const float T = 310.15; /* (K) */
const float INF = 10000;
/* Multiloop constants */
const float MC = 10.1; /* Multiloop initiation constant */
const float MB = -0.3; /* Multiloop bases don't paired constant */
const float MI = -0.3; /* Multiloop branchs of helice constant */

/* Variables Delaration */
string RNASeq;
float t;

/* Energy Matrix Declaration*/
float **W; // Matrix of total energy on RNA structure
float **V; // Matrix of energy for case of i, j basepair region
float **StackMatrix;
float **HairpinMatrix;
float **BulgeMatrix;
float **InternalMatrix;
float **MultiMatrix;
float ***MultiMatrix_1;
float ***MultiMatrix_2;

/* String Traceback Matrix Declaration */
string **W_TB;
string **V_TB;
string **Stack_TB;
string **Hairpin_TB;
string **Bulge_TB;
string **Internal_TB;
string **Multi_TB;
string ***Multi_1_TB;
string **Multi_2_TB;

/* Score for basepair */
map <string, float> BPScore;

/* Energy Maps Declaration */
map <string, float> Helices;
map <string, float> TermiMis;
map <string, float> InterTermiMis;
map <string, float> Inter2x3TermiMis;
map <string, float> Danglings;
// Hairpin loop special cases
map <string, float> HairpinSpecial;
// Internal loop special cases
map <string, float> Inter1x1;
map <string, float> Inter2x1;
map <string, float> Inter2x2;
// Loops Matrix: [size][]: [0]-internal, [1]-bulge; [2]-hairpin
float Loops[31][3];

/* Functions Declaration */
// Convert string to all uppercase

// Reading energy parameters of each basepair from file
void CallingPars();
void GettingBasepairPars(string FileName, map <string, float>& MapName);
void GettingLoopsPars(string FileName, float Array[][3]);

// Getting special cases of Internal loop (1x1, 2x1, 2x2)
inline string Inter_Model(int i, int j, int p, int q, int m, int x, int n, int y);

// Filling Matrix
float **MakingEnergyMatrix(float **Matrix, int size);
string **MakingStringMatrix(string **Matrix, int size);
float ***Making3DMatrix(float ***a, int size); // For MultiMatrix_1 and 2
string ***MakingString3DMatrix(string ***a, int size); // For MultiMatrix_TB_1
void GettingMatrix(string RNASeq);
float CalculatingMFE(string RNASeq);
float FillingV_Matrix(int i, int j);
float FillingHairpinMatrix(int i, int j);
float FillingStackMatrix(int i, int j);
float FillingBulgeMatrix(int i, int j);
float FillingInternalMatrix(int i, int j);
string Inter_Model(int i, int j, int p, int q, int m, int x, int n, int y);
float FillingMultiMatrix(int i, int j);
float FillingMulti_1_Matrix(int p, int q, int j);
float *FillingMulti_2_Matrix(int i, int p);
inline float FreeEnergy(float H, float E);
// Determining RNA structure
void RNAfolding(string RNASeq);

/**********************************************************************************/
/* MAIN PROGRAM */
/**********************************************************************************/
void main(){

	// Input RNA sequence
	while (RNASeq != "exit"){

		cout << "Input RNA sequence: " << endl;
		cin >> RNASeq;

		if (RNASeq == "exit"){
			exit(0);
		}
		else{
			transform(RNASeq.begin(), RNASeq.end(), RNASeq.begin(), (int(*)(int))toupper);
		}

		cout << endl;

		// Calculating RNA structure and its energy
		RNAfolding(RNASeq);
		cout << endl;
		cout << "*************************************************************" << endl;
	}
}

/**********************************************************************************/
/* FUNCTIONS */
/**********************************************************************************/
/* Convert string to all uppercase */

/* Reading energy parameters of each basepair from file*/
void GettingBasepairPars(string FileName, map <string, float>& MapName){
	string type;
	string arrow;
	float energy;

	// Opening energy file
	fstream inFile;
	inFile.open(FileName);

	// Check for Error opening File
	if (inFile.fail()){
		cerr << "Error Opening File " << FileName << endl;
		exit(1);
	}

	// Getting energy from type of sequence
	while (!inFile.eof()){
		inFile >> type >> arrow >> energy;
		MapName.insert(map<string, float>::value_type(type, energy)); // Insert parameters to Energy Map
	}
}

/* Reading loop energy for size from file*/
void GettingLoopsPars(string FileName, float Array[][3]){
	string size;
	string internalloop;
	string bulgeloop;
	string hairpinloop;
	int i = 0;

	// Opening energy file
	fstream inFile;
	inFile.open(FileName);

	// Check for Error opening File
	if (inFile.fail()){
		cerr << "Error Opening File " << FileName << endl;
		exit(1);
	}

	// Getting energy from loop size
	while (!inFile.eof()){
		inFile >> size >> internalloop >> bulgeloop >> hairpinloop;
		Array[i][0] = float(atof(internalloop.c_str()));
		Array[i][1] = float(atof(bulgeloop.c_str()));
		Array[i][2] = float(atof(hairpinloop.c_str()));
		i++;
	}
}

/* Calling all parameters from files */
void CallingPars(){
	GettingBasepairPars("BasepairScore.txt", BPScore);

	// Gibbs energy
	GettingBasepairPars("Watson Crick helices.txt", Helices);
	GettingBasepairPars("Terminal mismatchs.txt", TermiMis);
	GettingBasepairPars("Internal_Terminal mismatchs.txt", InterTermiMis);
	GettingBasepairPars("Internal_2x3_Terminal mismatchs.txt", Inter2x3TermiMis);
	GettingBasepairPars("Dangling_ends.txt", Danglings);
	GettingBasepairPars("HairpinSpecial.txt", HairpinSpecial);
	GettingBasepairPars("inter-1x1.txt", Inter1x1);
	GettingBasepairPars("inter-2x1.txt", Inter2x1);
	GettingBasepairPars("inter-2x2.txt", Inter2x2);
	GettingLoopsPars("Loops.txt", Loops);
}

/* Getting Matrixs */
float **MakingEnergyMatrix(float **Matrix, int size){
	Matrix = new float*[size];
	for (int i = 0; i < size; i++){
		Matrix[i] = new float[size];
	}

	return Matrix;
}

float ***Making3DMatrix(float ***a, int size){
	a = new float**[size];

	for (int i = 0; i < size; i++){
		a[i] = new float*[size];
	}

	for (int i = 0; i < size; i++){
		for (int j(0); j < size; j++){
			a[i][j] = new float[size];
		}
	}

	for (int i(0); i < size; i++){
		for (int j(0); j < size; j++){
			for (int k(0); k < size; k++){
				a[i][j][k] = NULL;
			}
		}
	}

	return a;
}

string **MakingStringMatrix(string **Matrix, int size){
	Matrix = new string*[size];
	for (int i = 0; i < size; i++){
		Matrix[i] = new string[size];
	}

	return Matrix;
}

string ***MakingString3DMatrix(string ***a, int size){
	a = new string**[size];

	for (int i = 0; i < size; i++){
		a[i] = new string*[size];
	}

	for (int i = 0; i < size; i++){
		for (int j(0); j < size; j++){
			a[i][j] = new string[size];
		}
	}

	return a;
}

void GettingMatrix(string RNAseq){
	int RNAlen = RNAseq.size();

	W = MakingEnergyMatrix(W, RNAlen);
	V = MakingEnergyMatrix(V, RNAlen);
	StackMatrix = MakingEnergyMatrix(StackMatrix, RNAlen);
	HairpinMatrix = MakingEnergyMatrix(HairpinMatrix, RNAlen);
	BulgeMatrix = MakingEnergyMatrix(BulgeMatrix, RNAlen);
	InternalMatrix = MakingEnergyMatrix(InternalMatrix, RNAlen);
	MultiMatrix = MakingEnergyMatrix(MultiMatrix, RNAlen);
	MultiMatrix_1 = Making3DMatrix(MultiMatrix_1, RNAlen);
	MultiMatrix_2 = Making3DMatrix(MultiMatrix_2, RNAlen);

	W_TB = MakingStringMatrix(W_TB, RNAlen);
	V_TB = MakingStringMatrix(V_TB, RNAlen);
	Stack_TB = MakingStringMatrix(Stack_TB, RNAlen);
	Hairpin_TB = MakingStringMatrix(Hairpin_TB, RNAlen);
	Bulge_TB = MakingStringMatrix(Bulge_TB, RNAlen);
	Internal_TB = MakingStringMatrix(Internal_TB, RNAlen);
	Multi_TB = MakingStringMatrix(Multi_TB, RNAlen);
	Multi_1_TB = MakingString3DMatrix(Multi_1_TB, RNAlen);
	Multi_2_TB = MakingStringMatrix(Multi_2_TB, RNAlen);
}

/* Determining structure of RNA may be folded */
void RNAfolding(string RNASeq){
	int RNAlen = RNASeq.size();

	// Calculating MFE and structure of RNA
	float MFE = CalculatingMFE(RNASeq);
	string structure = W_TB[0][RNAlen - 1];

	cout << "RNA structure: " << endl;
	cout << RNASeq << endl;
	cout << structure << endl;
	cout << "RNA length: " << RNAlen << endl;
	cout << "MFE = " << setprecision(4) << MFE << endl;

	delete[] V, V_TB, W, W_TB, StackMatrix, HairpinMatrix, BulgeMatrix, InternalMatrix, MultiMatrix, MultiMatrix_1, MultiMatrix_2,
			Stack_TB, Hairpin_TB, Bulge_TB, Internal_TB, Multi_TB, Multi_1_TB, Multi_2_TB;;
}
/* Calculating the minimum energy of RNA strucrture may abe folded */
float CalculatingMFE(string RNASeq){
	int RNAlen = RNASeq.size();
	int i = 0;

	// Time begin the program
	float start = clock();

	// Calling Pars From File
	CallingPars();

	// Creating matrix spaces
	GettingMatrix(RNASeq);

	// Filling V matrix
	for (int i = RNAlen - 5; i >= 0; --i){
		for (int j = i + 4; j < RNAlen; ++j){
			if (BPScore[RNASeq.substr(i, 1) + RNASeq.substr(j, 1)] == 1){
				// Calculating V(i,j) value
				V[i][j] = FillingV_Matrix(i, j);
			}
			else{
				V[i][j] = INF;
			}
		}
	}

	// Filling W matricx (Total energy matrix)
	for (int j = 0; j <= 3; ++j){
		W[i][j] = INF;
		string dots(j - i + 1, '.');
		W_TB[i][j] = dots;
	}

	for (int j = i + 4; j < RNAlen; ++j){
		// Case j not paired
		W[i][j] = W[i][j - 1];
		W_TB[i][j] = W_TB[i][j - 1] + ".";

		// Case i,j paired
		if (BPScore[RNASeq.substr(i, 1) + RNASeq.substr(j, 1)] == 1){
			if (W[i][j] > V[i][j]){
				W[i][j] = V[i][j];
				W_TB[i][j] = V_TB[i][j];
			}
		}

		// Case j paired with k between (i,j)
		for (int k = i + 1; k <= j - 4; ++k){
			if (BPScore[RNASeq.substr(k, 1) + RNASeq.substr(j, 1)] == 1){
				float W_min = W[i][k - 1] + V[k][j];
				if (k > 0){
					W_min += Danglings[RNASeq.substr(k - 1, 1) + "-" + RNASeq.substr(j, 1) + RNASeq.substr(k, 1)];
				}
				if (j < (RNAlen - 1)){
					W_min += Danglings[RNASeq.substr(j, 1) + RNASeq.substr(k, 1) + "-" + RNASeq.substr(j + 1, 1)];
				}
				if (W[i][j] > W_min){
					W[i][j] = W_min;
					W_TB[i][j] = W_TB[i][k - 1] + V_TB[k][j];
				}
			}
		}
	}

	// Time for finish the program
	float finish = clock();

	float TimeTotal = (float)(finish - start) / CLOCKS_PER_SEC;

	cout << "Time: " << setprecision(3) << TimeTotal << " (s)" << endl;

	return W[0][RNAlen - 1];
}

float FillingV_Matrix(int i, int j){
	float Vij = INF;
	float Stack = INF;
	float Bulge = INF;
	float Internal = INF;
	float Multi = INF;

	// Hairpin
	float Hairpin = FillingHairpinMatrix(i, j);

	// Stack
	if (j - i >= 6){
		Stack = FillingStackMatrix(i, j);
	}
	//Bulgeloop
	if (j - i >= 7){
		Bulge = FillingBulgeMatrix(i, j);
	}
	//Internal loop
	if (j - i >= 8){
		Internal = FillingInternalMatrix(i, j);
	}

	// Multiloop
	if (j - i >= 11){
		Multi = FillingMultiMatrix(i, j);
	}

	float EnergyGroup[5] = { Stack, Hairpin, Bulge, Internal, Multi };
	string StringGroup[5] = { Stack_TB[i][j], Hairpin_TB[i][j], Bulge_TB[i][j], Internal_TB[i][j], Multi_TB[i][j] };

	for (int k = 0; k < 5; k++){
		if (EnergyGroup[k] < Vij){
			Vij = EnergyGroup[k];
			V_TB[i][j] = StringGroup[k];
		}
	}

	return Vij;
}

/* Determining type of structure of RNA may be folded */
float FillingHairpinMatrix(int i, int j){
	float HairpinE;

	string dots(j - i - 1, '.');
	Hairpin_TB[i][j] = "(" + dots + ")";
	int RNAlen = RNASeq.size();

	// Calculating for special hairpin loop
	if (j - i - 1 == 3 || j - i - 1 == 4 || j - i - 1 == 6){ // Case for 4 or 6 nu hairpinloop
		for (map<string, float>::iterator it = HairpinSpecial.begin(); it != HairpinSpecial.end(); it++){
			string key = it->first;
			if (RNASeq.substr(i, j - i + 1) == key){
				HairpinE = HairpinSpecial[key];
				return HairpinE;
			}
		}
	}

	// Other hairpin loops
	// Calculating all parameters
	float Termi_mismatch;
	float Mismatch_E_bonus = 0;
	float GU_close_bonus = 0; // Case that basepair(i, j) is GU
	float C_loops_E_bonus = 0; // Case loop have nu C
	float AU_end_penalty = 0; // Case basepair(i, j) is AU

	// C_loop_E_bonus: loop with all nucleotide C
	if (RNASeq.substr(i + 1, j - i - 1).find("A") == string::npos && RNASeq.substr(i + 1, j - i - 1).find("G") == string::npos
		&& RNASeq.substr(i + 1, j - i - 1).find("U") == string::npos){
		if ((j - i) == 4){
			C_loops_E_bonus += 1.5;
		}
		else{
			C_loops_E_bonus += 0.3*(j - i - 1) + 1.6;
		}
	}

	// Case if hairpin have 3 nu in loop
	if (j - i == 4){ // Case for 3 nu hairpinloop
		HairpinE = Loops[3][2] + C_loops_E_bonus;
		return HairpinE;
	}

	// Terminal mismatch
	string subij = RNASeq.substr(i, 1) + RNASeq.substr(j, 1);
	string subi1j1 = RNASeq.substr(i + 1, 1) + RNASeq.substr(j - 1, 1);
	Termi_mismatch = TermiMis[subij + "-" + subi1j1];

	// Mismatch Energy bonus
	if (subi1j1 == "UU" || subi1j1 == "GA"){
		Mismatch_E_bonus += -0.9;
	}
	else if (subi1j1 == "GG"){
		Mismatch_E_bonus += -0.8;
	}

	// AU_end_penalty
	if (subij == "AU" || subij == "UA"){
		AU_end_penalty += 0.45;
	}
	// GU_close_bonus
	else if ((i >= 2 && subij == "GU" && RNASeq.substr(i - 2, 2) == "GG")
		|| (j < RNAlen - 3 && subij == "UG" && RNASeq.substr(j + 1, 2) == "GG")){
		GU_close_bonus += -2.2;
	}

	// Calculating hairpin loop energy
	HairpinE = Termi_mismatch + Mismatch_E_bonus + GU_close_bonus + C_loops_E_bonus + AU_end_penalty;
	if ((j - i) <= 31){
		HairpinE += Loops[j - i - 1][2];
	}
	else{
		HairpinE += Loops[9][2] + 1.75*R*T*log((j - i - 1) / 9);
	}

	return HairpinE;
}

float FillingStackMatrix(int i, int j){
	float StackE = INF;
	float HeliceE;

	if (BPScore[RNASeq.substr(i + 1, 1) + RNASeq.substr(j - 1, 1)] == 1){

		if (RNASeq.substr(i, 4) == "GGUC" && RNASeq.substr(j - 3, 4) == "GGUC"){
			HeliceE = -4.12;
			StackE = HeliceE + V[i + 3][j - 3];
			Stack_TB[i][j] = "(((" + V_TB[i + 3][j - 3] + ")))";
			return StackE;
		}
		else{
			HeliceE = Helices[RNASeq.substr(i, 1) + RNASeq.substr(j, 1) + "-" + RNASeq.substr(i + 1, 1) + RNASeq.substr(j - 1, 1)];
			StackE = HeliceE + V[i + 1][j - 1];
			Stack_TB[i][j] = "(" + V_TB[i + 1][j - 1] + ")";
			return StackE;
		}
	}

	return StackE;
}

float FillingBulgeMatrix(int i, int j){
	float BulgeE = INF;
	float BulgeEmin = INF;
	string structure(j - i - 1, '.');

	// 2 Cases: 1. Bulgeloop in 3'end
	//			2. Bulgeloop in 5'end

	// Case 1:
	for (int q = i + 5; q <= j - 2; ++q){
		if (BPScore[RNASeq.substr(i + 1, 1) + RNASeq.substr(q, 1)] == 1){
			// Calculating Bulgeloop Energy
			if (j - q == 2){ // Case loop have only 1 nu
				float Special_C = 0;
				if ((q <= j - 3 && RNASeq.substr(q + 1, 3) == "CCC") || RNASeq.substr(q - 1, 3) == "CCC"){
					Special_C = -0.9;
				}

				BulgeE = Loops[1][1] + Special_C
					+ Helices[RNASeq.substr(i, 1) + RNASeq.substr(j, 1) + "-" + RNASeq.substr(i + 1, 1) + RNASeq.substr(q, 1)]
					+ V[i + 1][q];
			}
			else{ // Case loop have more than 1 nu
				float AU_end_penalty = 0; // Case of bulge loop closed by AU
				float GU_end_penalty = 0; // Case of bulge loop closed by GU
				string subij = RNASeq.substr(i, 1) + RNASeq.substr(j, 1);
				string subqi1 = RNASeq.substr(q, 1) + RNASeq.substr(i + 1, 1);

				if (subij == "AU" || subij == "UA"){
					AU_end_penalty += 0.45;
				}
				else if (subij == "GU" || subij == "UG"){
					GU_end_penalty += 0.45;
				}

				if (subqi1 == "AU" || subqi1 == "UA"){
					AU_end_penalty += 0.45;
				}
				else if (subqi1 == "GU" || subqi1 == "UG"){
					GU_end_penalty += 0.45;
				}

				BulgeE = AU_end_penalty + GU_end_penalty + V[i + 1][q];
				if (j - q <= 31){
					BulgeE += Loops[j - q - 1][1];
				}
				else{
					BulgeE += Loops[6][1] + 1.75*R*T*log((j - q - 1) / 6);
				}
			}

			if (BulgeEmin > BulgeE){
				BulgeEmin = BulgeE;
				string dots(j - q - 1, '.');
				structure = "(" + V_TB[i + 1][q] + dots + ")";
			}
		}
	}

	// Case 2:
	for (int p = i + 2; p <= j - 5; ++p){
		if (BPScore[RNASeq.substr(p, 1) + RNASeq.substr(j - 1, 1)] == 1){
			// Calculating Bulgeloop Energy
			if (p - i == 2){ // Case loop have only 1 nu
				float Special_C = 0;
				if ((p >= 3 && RNASeq.substr(p - 3, 3) == "CCC") || RNASeq.substr(p - 1, 3) == "CCC"){
					Special_C = -0.9;
				}
				BulgeE = Loops[1][1] + Special_C
					+ Helices[RNASeq.substr(i, 1) + RNASeq.substr(j, 1) + "-" + RNASeq.substr(p, 1) + RNASeq.substr(j - 1, 1)]
					+ V[p][j - 1];
			}
			else{
				float AU_end_penalty = 0;
				float GU_end_penalty = 0;
				string subij = RNASeq.substr(i, 1) + RNASeq.substr(j, 1);
				string subj1p = RNASeq.substr(j - 1, 1) + RNASeq.substr(p, 1);
				if (subij == "AU" || subij == "UA"){
					AU_end_penalty += 0.5;
				}
				else if (subij == "GU" || subij == "UG"){
					GU_end_penalty += 0.47;
				}

				if (subj1p == "AU" || subj1p == "UA"){
					AU_end_penalty += 0.5;
				}
				else if (subj1p == "GU" || subj1p == "UG"){
					GU_end_penalty += 0.47;
				}

				BulgeE = AU_end_penalty + GU_end_penalty + V[p][j - 1];
				if (p - i <= 31){
					BulgeE += Loops[p - i - 1][1];
				}
				else{
					BulgeE += Loops[6][1] + 1.75*R*T*log((p - i - 1) / 6);
				}
			}

			if (BulgeEmin > BulgeE){
				BulgeEmin = BulgeE;
				string dots(p - i - 1, '.');
				structure = "(" + dots + V_TB[p][j - 1] + ")";
			}
		}
	}

	Bulge_TB[i][j] = structure;
	return BulgeEmin;
}

float FillingInternalMatrix(int i, int j){
	float InternalE;
	float InternalEMin = INF;
	float AU_GU_closure;
	float TerMismatch_1; // Mismatch (i,j) end
	float TerMismatch_2; // Mismatch (p,q) end
	string structure(j - i - 1, '.');

	for (int p = i + 2; p < j - 6; ++p){
		for (int q = p + 4; q <= j - 2; ++q){
			if (BPScore[RNASeq.substr(p, 1) + RNASeq.substr(q, 1)] == 1){

				// Calculating energy of Internal loop

				//Case 1: Loop 1x1
				if ((p - i) == 2 && (j - q) == 2){
					InternalE = Inter1x1[Inter_Model(i, j, p, q, p - 1, 1, q + 1, 1)] + V[p][q];
				}

				// Case 2: Loop 2x1
				else if ((p - i) == 2 && (j - q) == 3){
					InternalE = Inter2x1[Inter_Model(i, j, p, q, p - 1, 1, q + 1, 2)] + V[p][q];
				}

				else if ((p - i) == 3 && (j - q) == 2){
					InternalE = Inter2x1[Inter_Model(q, p, j, i, q + 1, 1, p - 2, 2)] + V[p][q];
				}

				// Case 3: Loop 2x2
				else if ((p - i) == 3 && (j - q) == 3){
					InternalE = Inter2x2[Inter_Model(i, j, p, q, p - 2, 2, q + 1, 2)] + V[p][q];
				}

				// Case 4: Others
				else{
					// Calculating parameters
					// AU_GU_closure
					AU_GU_closure = 0;
					string subij = RNASeq.substr(i, 1) + RNASeq.substr(j, 1);
					string subpq = RNASeq.substr(q, 1) + RNASeq.substr(p, 1);
					if (subij == "AU" || subij == "GU"
						|| subij == "UA" || subij == "UG"){
						AU_GU_closure += 0.7;
					}
					if (subpq == "AU" || subpq == "GU"
						|| subpq == "UA" || subpq == "UG"){
						AU_GU_closure += 0.7;
					}

					// TerMismatch_1 and TerMismatch_2
					if ((p - i - 1) == 1 || (j - q - 1) == 1){
						TerMismatch_1 = 0;
						TerMismatch_2 = 0;
					}
					else if (((p - i) == 3 && (j - q) == 4) || ((p - i - 1) == 3 && (j - q - 1) == 2)){
						TerMismatch_1 = Inter2x3TermiMis[subij + "-" + RNASeq.substr(i + 1, 1) + RNASeq.substr(j - 1, 1)];
						TerMismatch_2 = Inter2x3TermiMis[subpq + "-" + RNASeq.substr(q + 1, 1) + RNASeq.substr(p - 1, 1)];
					}
					else{
						TerMismatch_1 = InterTermiMis[subij + "-" + RNASeq.substr(i + 1, 1) + RNASeq.substr(j - 1, 1)];
						TerMismatch_2 = InterTermiMis[subpq + "-" + RNASeq.substr(q + 1, 1) + RNASeq.substr(p - 1, 1)];
					}

					// Calculating Internal loop energy
					InternalE = 0.6 * abs(p - i - j + q) + TerMismatch_1 + TerMismatch_2 + AU_GU_closure + V[p][q];

					if ((p - i + j - q) >= 6 && (p - i + j - q) <= 32){
						InternalE += Loops[p - i + j - q - 2][0];
					}
					else{
						InternalE += Loops[6][0] + 1.08 * log((p - i + j - q - 2) / 6);
					}
				}

				if (InternalEMin > InternalE){
					InternalEMin = InternalE;
					string dots_1(p - i - 1, '.');
					string dots_2(j - q - 1, '.');
					structure = "(" + dots_1 + V_TB[p][q] + dots_2 + ")";
				}
			}
		}
	}

	Internal_TB[i][j] = structure;
	return InternalEMin;
}

string Inter_Model(int i, int j, int p, int q, int m, int x, int n, int y){
	// m: Position of nucleotide don't pair on the top sequence
	// n: Position of nucleotide don't pair on the bottom sequence
	return RNASeq.substr(i, 1) + RNASeq.substr(j, 1) + "-" + RNASeq.substr(m, x) + "-" + RNASeq.substr(n, y) + "-" + RNASeq.substr(p, 1) + RNASeq.substr(q, 1);
}

float FillingMultiMatrix(int i, int j){
	float MultiE;
	float MultiEMin = INF;
	string structure(j - i - 1, '.');

	for (int p = i + 6; p <= j - 5; ++p){
		for (int q = p + 4; q <= j - 1; ++q){
			if (BPScore[RNASeq.substr(p, 1) + RNASeq.substr(q, 1)] == 1){
				//MultiE = Initiation Energy + Stack Energy

				if (MultiMatrix_1[p][q][j] == NULL){
					MultiMatrix_1[p][q][j] = FillingMulti_1_Matrix(p, q, j);
				}
				if (MultiMatrix_2[i][p][0] == NULL){
					float *Pars = new float[3];
					Pars = FillingMulti_2_Matrix(i, p);
					MultiMatrix_2[i][p][0] = Pars[0];
					MultiMatrix_2[i][p][1] = Pars[1];
					MultiMatrix_2[i][p][2] = Pars[2];
					delete[] Pars;
				}

				MultiE = MC + MultiMatrix_1[p][q][j] + MultiMatrix_2[i][p][0];

				if (j - q > 2){
					MultiE += Danglings[RNASeq.substr(j - 1, 1) + "-" + RNASeq.substr(i, 1) + RNASeq.substr(j, 1)];
				}
				if (p - MultiMatrix_2[i][p][1] > 2){
					MultiE += Danglings[RNASeq.substr(p - 1, 1) + "-" + RNASeq.substr(q, 1) + RNASeq.substr(p, 1)];
				}
				if (MultiMatrix_2[i][p][2] - i > 1){
					MultiE += Danglings[RNASeq.substr(i, 1) + RNASeq.substr(j, 1) + "-" + RNASeq.substr(i + 1, 1)];
				}

				if (MultiEMin > MultiE){
					MultiEMin = MultiE;
					structure = "(" + Multi_2_TB[i][p] + Multi_1_TB[p][q][j] + ")";
				}
			}
		}
	}

	Multi_TB[i][j] = structure;
	return MultiEMin;
}

float FillingMulti_1_Matrix(int p, int q, int j){
	string dots(j - q - 1, '.');
	Multi_1_TB[p][q][j] = V_TB[p][q] + dots;

	if (j - q - 1 >= 1){
		return V[p][q] + MB*(j - q - 1) + MI + Danglings[RNASeq.substr(q, 1) + RNASeq.substr(p, 1) + "-" + RNASeq.substr(q + 1, 1)];
	}
	else{
		return V[p][q] + MB*(j - q - 1) + MI;
	}
}

float *FillingMulti_2_Matrix(int i, int p){
	float ML_2_1 = INF;
	float ML_2_m = INF;
	float *Pars = new float[3];
	Pars[0] = INF;
	Pars[1] = -1;
	Pars[2] = -1;
	string dots(p - i - 1, '.');
	Multi_2_TB[i][p] = dots;

	for (int k = i + 1; k <= p - 5; k++){
		for (int n = k + 4; n <= p - 1; n++){
			if (BPScore[RNASeq.substr(k, 1) + RNASeq.substr(n, 1)] == 1){
				if (MultiMatrix_1[k][n][p] == NULL){
					MultiMatrix_1[k][n][p] = FillingMulti_1_Matrix(k, n, p);
				}
				if (MultiMatrix_2[i][k][0] == NULL){
					float *APars = new float[3];
					APars = FillingMulti_2_Matrix(i, k);
					MultiMatrix_2[i][k][0] = APars[0];
					MultiMatrix_2[i][k][1] = APars[1];
					MultiMatrix_2[i][k][2] = APars[2];
					delete[] APars;
				}

				ML_2_1 = MultiMatrix_1[k][n][p] + MB*(k - i - 1) + Danglings[RNASeq.substr(k - 1, 1) + "-" + RNASeq.substr(n, 1) + RNASeq.substr(k, 1)];
				ML_2_m = MultiMatrix_1[k][n][p] + MultiMatrix_2[i][k][0];

				if (k - MultiMatrix_2[i][k][1] > 2){
					ML_2_m += Danglings[RNASeq.substr(k - 1, 1) + "-" + RNASeq.substr(n, 1) + RNASeq.substr(k, 1)];
				}

				if (Pars[0] > ML_2_1){
					Pars[0] = ML_2_1;
					Pars[1] = n;
					Pars[2] = k;
				}

				if (Pars[0] > ML_2_m){
					Pars[0] = ML_2_m;
					Pars[1] = n;
					Pars[2] = MultiMatrix_2[i][k][2];
				}

				// Deteming structure of Multiloop
				if (Pars[0] == ML_2_1){
					string dots(k - i - 1, '.');
					Multi_2_TB[i][p] = dots + Multi_1_TB[k][n][p];
				}
				if (Pars[0] == ML_2_m){
					Multi_2_TB[i][p] = Multi_2_TB[i][k] + Multi_1_TB[k][n][p];
				}
			}
		}
	}

	return Pars;
}