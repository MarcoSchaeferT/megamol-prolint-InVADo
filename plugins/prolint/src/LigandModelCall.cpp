/*
 * LigandModelCall.cpp
 *
 * Copyright (C) 2010 by University of Stuttgart (VISUS).
 * All rights reserved.
 */

#include "stdafx.h"
#include "LigandModelCall.h"
#include "vislib/IllegalParamException.h"
#include "vislib/IllegalStateException.h"
#include "vislib/math/mathfunctions.h"
#include "vislib/OutOfRangeException.h"

using namespace megamol;
using namespace megamol::prolint;

/*
 * LigandModelCall::CallForGetData
 */
const unsigned int LigandModelCall::CallForGetData = 0;

/*
 * LigandModelCall::CallForGetExtent
 */
const unsigned int LigandModelCall::CallForGetExtent = 1;

/*
 * LigandModelCall::CallForGetDataSilent
 */
const unsigned int LigandModelCall::CallForGetDataSilent = 2;

/*
 * LigandModelCall::CallForGetClusterAndCentroidData
 */
const unsigned int LigandModelCall::CallForGetClusterAndCentroidData = 3;

/*
 * LigandModelCall::CallForCallForSetLMCStatus
 */
const unsigned int LigandModelCall::CallForSetLMCStatus = 4;

/*
 * LigandModelCall::CallForSetWebData
 */
const unsigned int LigandModelCall::CallForSetWebData = 5;


/*
 * LigandModelCall::LigandModelCall
 */
LigandModelCall::LigandModelCall(void)
    : AbstractGetDataCall(),
	currentLigand(0),
	currentModel(0),
	ligandCount(0),
	modelCount(NULL),
	currentGlobalMdlID(-1),
	globalMdlID_setByWeb(-1){

	this->functionalGroupsIDsToWord.resize(205);
    this->functionalGroupsIDsToWord[1] = "cation";
    this->functionalGroupsIDsToWord[2] = "anion";
    this->functionalGroupsIDsToWord[3] = "carbonyl compound";
    this->functionalGroupsIDsToWord[4] = "aldehyde";
    this->functionalGroupsIDsToWord[5] = "ketone";
    this->functionalGroupsIDsToWord[6] = "thiocarbonyl compound";
    this->functionalGroupsIDsToWord[7] = "thioaldehyde";
    this->functionalGroupsIDsToWord[8] = "thioketone";
    this->functionalGroupsIDsToWord[9] = "imine";
    this->functionalGroupsIDsToWord[10] = "hydrazone";
    this->functionalGroupsIDsToWord[11] = "semicarbazone";
    this->functionalGroupsIDsToWord[12] = "thiosemicarbazone";
    this->functionalGroupsIDsToWord[13] = "oxime";
    this->functionalGroupsIDsToWord[14] = "oxime ether";
    this->functionalGroupsIDsToWord[15] = "ketene";
    this->functionalGroupsIDsToWord[16] = "ketene acetal or derivative";
    this->functionalGroupsIDsToWord[17] = "carbonyl hydrate";
    this->functionalGroupsIDsToWord[18] = "hemiacetal";
    this->functionalGroupsIDsToWord[19] = "acetal";
    this->functionalGroupsIDsToWord[20] = "hemiaminal";
    this->functionalGroupsIDsToWord[21] = "aminal";
    this->functionalGroupsIDsToWord[22] = "hemithioaminal";
    this->functionalGroupsIDsToWord[23] = "thioacetal";
    this->functionalGroupsIDsToWord[24] = "enamine";
    this->functionalGroupsIDsToWord[25] = "enol";
    this->functionalGroupsIDsToWord[26] = "enol ether";
    this->functionalGroupsIDsToWord[27] = "hydroxy compound";
    this->functionalGroupsIDsToWord[28] = "alcohol";
    this->functionalGroupsIDsToWord[29] = "primary alcohol";
    this->functionalGroupsIDsToWord[30] = "secondary alcohol";
    this->functionalGroupsIDsToWord[31] = "tertiary alcohol";
    this->functionalGroupsIDsToWord[32] = "1,2-diol";
    this->functionalGroupsIDsToWord[33] = "1,2-aminoalcohol";
    this->functionalGroupsIDsToWord[34] = "phenol or hydroxyhetarene";
    this->functionalGroupsIDsToWord[35] = "1,2-diphenol";
    this->functionalGroupsIDsToWord[36] = "enediol";
    this->functionalGroupsIDsToWord[37] = "ether";
    this->functionalGroupsIDsToWord[38] = "dialkyl ether";
    this->functionalGroupsIDsToWord[39] = "alkyl aryl ether";
    this->functionalGroupsIDsToWord[40] = "diaryl ether";
    this->functionalGroupsIDsToWord[41] = "thioether";
    this->functionalGroupsIDsToWord[42] = "disulfide";
    this->functionalGroupsIDsToWord[43] = "peroxide";
    this->functionalGroupsIDsToWord[44] = "hydroperoxide";
    this->functionalGroupsIDsToWord[45] = "hydrazine derivative";
    this->functionalGroupsIDsToWord[46] = "hydroxylamine";
    this->functionalGroupsIDsToWord[47] = "amine";
    this->functionalGroupsIDsToWord[48] = "primary amine";
    this->functionalGroupsIDsToWord[49] = "primary aliphatic amine (alkylamine)";
    this->functionalGroupsIDsToWord[50] = "primary aromatic amine";
    this->functionalGroupsIDsToWord[51] = "secondary amine";
    this->functionalGroupsIDsToWord[52] = "secondary aliphatic amine (dialkylamine)";
    this->functionalGroupsIDsToWord[53] = "secondary aliphatic/aromatic amine (alkylarylamine)";
    this->functionalGroupsIDsToWord[54] = "secondary aromatic amine (diarylamine)";
    this->functionalGroupsIDsToWord[55] = "tertiary amine";
    this->functionalGroupsIDsToWord[56] = "tertiary aliphatic amine (trialkylamine)";
    this->functionalGroupsIDsToWord[57] = "tertiary aliphatic/aromatic amine (alkylarylamine)";
    this->functionalGroupsIDsToWord[58] = "tertiary aromatic amine (triarylamine)";
    this->functionalGroupsIDsToWord[59] = "quaternary ammonium salt";
    this->functionalGroupsIDsToWord[60] = "N-oxide";
    this->functionalGroupsIDsToWord[61] = "halogen derivative";
    this->functionalGroupsIDsToWord[62] = "alkyl halide";
    this->functionalGroupsIDsToWord[63] = "alkyl fluoride";
    this->functionalGroupsIDsToWord[64] = "alkyl chloride";
    this->functionalGroupsIDsToWord[65] = "alkyl bromide";
    this->functionalGroupsIDsToWord[66] = "alkyl iodide";
    this->functionalGroupsIDsToWord[67] = "aryl halide";
    this->functionalGroupsIDsToWord[68] = "aryl fluoride";
    this->functionalGroupsIDsToWord[69] = "aryl chloride";
    this->functionalGroupsIDsToWord[70] = "aryl bromide";
    this->functionalGroupsIDsToWord[71] = "aryl iodide";
    this->functionalGroupsIDsToWord[72] = "organometallic compound";
    this->functionalGroupsIDsToWord[73] = "organolithium compound";
    this->functionalGroupsIDsToWord[74] = "organomagnesium compound";
    this->functionalGroupsIDsToWord[75] = "carboxylic acid derivative";
    this->functionalGroupsIDsToWord[76] = "carboxylic acid";
    this->functionalGroupsIDsToWord[77] = "carboxylic acid salt";
    this->functionalGroupsIDsToWord[78] = "carboxylic acid ester";
    this->functionalGroupsIDsToWord[79] = "lactone";
    this->functionalGroupsIDsToWord[80] = "carboxylic acid amide";
    this->functionalGroupsIDsToWord[81] = "primary carboxylic acid amide";
    this->functionalGroupsIDsToWord[82] = "secondary carboxylic acid amide";
    this->functionalGroupsIDsToWord[83] = "tertiary carboxylic acid amide";
    this->functionalGroupsIDsToWord[84] = "lactam";
    this->functionalGroupsIDsToWord[85] = "carboxylic acid hydrazide";
    this->functionalGroupsIDsToWord[86] = "carboxylic acid azide";
    this->functionalGroupsIDsToWord[87] = "hydroxamic acid";
    this->functionalGroupsIDsToWord[88] = "carboxylic acid amidine";
    this->functionalGroupsIDsToWord[89] = "carboxylic acid amidrazone";
    this->functionalGroupsIDsToWord[90] = "carbonitrile";
    this->functionalGroupsIDsToWord[91] = "acyl halide";
    this->functionalGroupsIDsToWord[92] = "acyl fluoride";
    this->functionalGroupsIDsToWord[93] = "acyl chloride";
    this->functionalGroupsIDsToWord[94] = "acyl bromide";
    this->functionalGroupsIDsToWord[95] = "acyl iodide";
    this->functionalGroupsIDsToWord[96] = "acyl cyanide";
    this->functionalGroupsIDsToWord[97] = "imido ester";
    this->functionalGroupsIDsToWord[98] = "imidoyl halide";
    this->functionalGroupsIDsToWord[99] = "thiocarboxylic acid derivative";
    this->functionalGroupsIDsToWord[100] = "thiocarboxylic acid";
    this->functionalGroupsIDsToWord[101] = "thiocarboxylic acid ester";
    this->functionalGroupsIDsToWord[102] = "thiolactone";
    this->functionalGroupsIDsToWord[103] = "thiocarboxylic acid amide";
    this->functionalGroupsIDsToWord[104] = "thiolactam";
    this->functionalGroupsIDsToWord[105] = "imidothioester";
    this->functionalGroupsIDsToWord[106] = "oxo(het)arene";
    this->functionalGroupsIDsToWord[107] = "thioxo(het)arene";
    this->functionalGroupsIDsToWord[108] = "imino(het)arene";
    this->functionalGroupsIDsToWord[109] = "orthocarboxylic acid derivative";
    this->functionalGroupsIDsToWord[110] = "orthoester";
    this->functionalGroupsIDsToWord[111] = "amide acetal";
    this->functionalGroupsIDsToWord[112] = "carboxylic acid anhydride";
    this->functionalGroupsIDsToWord[113] = "carboxylic acid imide";
    this->functionalGroupsIDsToWord[114] = "carboxylic acid imide, N-unsubstituted";
    this->functionalGroupsIDsToWord[115] = "carboxylic acid imide, N-substituted";
    this->functionalGroupsIDsToWord[116] = "CO2 derivative (general)";
    this->functionalGroupsIDsToWord[117] = "carbonic acid derivative";
    this->functionalGroupsIDsToWord[118] = "carbonic acid monoester";
    this->functionalGroupsIDsToWord[119] = "carbonic acid diester";
    this->functionalGroupsIDsToWord[120] = "carbonic acid ester halide (alkyl/aryl haloformate)";
    this->functionalGroupsIDsToWord[121] = "thiocarbonic acid derivative";
    this->functionalGroupsIDsToWord[122] = "thiocarbonic acid monoester";
    this->functionalGroupsIDsToWord[123] = "thiocarbonic acid diester";
    this->functionalGroupsIDsToWord[124] = "thiocarbonic acid ester halide (alkyl/aryl halothioformate";
    this->functionalGroupsIDsToWord[125] = "carbamic acid derivative";
    this->functionalGroupsIDsToWord[126] = "carbamic acid";
    this->functionalGroupsIDsToWord[127] = "carbamic acid ester (urethane)";
    this->functionalGroupsIDsToWord[128] = "carbamic acid halide (haloformic acid amide)";
    this->functionalGroupsIDsToWord[129] = "thiocarbamic acid derivative";
    this->functionalGroupsIDsToWord[130] = "thiocarbamic acid";
    this->functionalGroupsIDsToWord[131] = "thiocarbamic acid ester";
    this->functionalGroupsIDsToWord[132] = "thiocarbamic acid halide (halothioformic acid amide)";
    this->functionalGroupsIDsToWord[133] = "urea";
    this->functionalGroupsIDsToWord[134] = "isourea";
    this->functionalGroupsIDsToWord[135] = "thiourea";
    this->functionalGroupsIDsToWord[136] = "isothiourea";
    this->functionalGroupsIDsToWord[137] = "guanidine";
    this->functionalGroupsIDsToWord[138] = "semicarbazide";
    this->functionalGroupsIDsToWord[139] = "thiosemicarbazide";
    this->functionalGroupsIDsToWord[140] = "azide";
    this->functionalGroupsIDsToWord[141] = "azo compound";
    this->functionalGroupsIDsToWord[142] = "diazonium salt";
    this->functionalGroupsIDsToWord[143] = "isonitrile";
    this->functionalGroupsIDsToWord[144] = "cyanate";
    this->functionalGroupsIDsToWord[145] = "isocyanate";
    this->functionalGroupsIDsToWord[146] = "thiocyanate";
    this->functionalGroupsIDsToWord[147] = "isothiocyanate";
    this->functionalGroupsIDsToWord[148] = "carbodiimide";
    this->functionalGroupsIDsToWord[149] = "nitroso compound";
    this->functionalGroupsIDsToWord[150] = "nitro compound";
    this->functionalGroupsIDsToWord[151] = "nitrite";
    this->functionalGroupsIDsToWord[152] = "nitrate";
    this->functionalGroupsIDsToWord[153] = "sulfuric acid derivative";
    this->functionalGroupsIDsToWord[154] = "sulfuric acid";
    this->functionalGroupsIDsToWord[155] = "sulfuric acid monoester";
    this->functionalGroupsIDsToWord[156] = "sulfuric acid diester";
    this->functionalGroupsIDsToWord[157] = "sulfuric acid amide ester";
    this->functionalGroupsIDsToWord[158] = "sulfuric acid amide";
    this->functionalGroupsIDsToWord[159] = "sulfuric acid diamide";
    this->functionalGroupsIDsToWord[160] = "sulfuryl halide";
    this->functionalGroupsIDsToWord[161] = "sulfonic acid derivative";
    this->functionalGroupsIDsToWord[162] = "sulfonic acid";
    this->functionalGroupsIDsToWord[163] = "sulfonic acid ester";
    this->functionalGroupsIDsToWord[164] = "sulfonamide";
    this->functionalGroupsIDsToWord[165] = "sulfonyl halide";
    this->functionalGroupsIDsToWord[166] = "sulfone";
    this->functionalGroupsIDsToWord[167] = "sulfoxide";
    this->functionalGroupsIDsToWord[168] = "sulfinic acid derivative";
    this->functionalGroupsIDsToWord[169] = "sulfinic acid";
    this->functionalGroupsIDsToWord[170] = "sulfinic acid ester";
    this->functionalGroupsIDsToWord[171] = "sulfinic acid halide";
    this->functionalGroupsIDsToWord[172] = "sulfinic acid amide";
    this->functionalGroupsIDsToWord[173] = "sulfenic acid derivative";
    this->functionalGroupsIDsToWord[174] = "sulfenic acid";
    this->functionalGroupsIDsToWord[175] = "sulfenic acid ester";
    this->functionalGroupsIDsToWord[176] = "sulfenic acid halide";
    this->functionalGroupsIDsToWord[177] = "sulfenic acid amide";
    this->functionalGroupsIDsToWord[178] = "thiol (sulfanyl compound)";
    this->functionalGroupsIDsToWord[179] = "alkylthiol";
    this->functionalGroupsIDsToWord[180] = "arylthiol (thiophenol)";
    this->functionalGroupsIDsToWord[181] = "phosphoric acid derivative";
    this->functionalGroupsIDsToWord[182] = "phosphoric acid";
    this->functionalGroupsIDsToWord[183] = "phosphoric acid ester";
    this->functionalGroupsIDsToWord[184] = "phosphoric acid halide";
    this->functionalGroupsIDsToWord[185] = "phosphoric acid amide";
    this->functionalGroupsIDsToWord[186] = "thiophosphoric acid derivative";
    this->functionalGroupsIDsToWord[187] = "thiophosphoric acid";
    this->functionalGroupsIDsToWord[188] = "thiophosphoric acid ester";
    this->functionalGroupsIDsToWord[189] = "thiophosphoric acid halide";
    this->functionalGroupsIDsToWord[190] = "thiophosphoric acid amide";
    this->functionalGroupsIDsToWord[191] = "phosphonic acid derivative";
    this->functionalGroupsIDsToWord[192] = "phosphonic acid";
    this->functionalGroupsIDsToWord[193] = "phosphonic acid ester";
    this->functionalGroupsIDsToWord[194] = "phosphine";
    this->functionalGroupsIDsToWord[195] = "phosphine oxide";
    this->functionalGroupsIDsToWord[196] = "boronic acid derivative";
    this->functionalGroupsIDsToWord[197] = "boronic acid";
    this->functionalGroupsIDsToWord[198] = "boronic acid ester";
    this->functionalGroupsIDsToWord[199] = "alkene";
    this->functionalGroupsIDsToWord[200] = "alkyne";
    this->functionalGroupsIDsToWord[201] = "aromatic compound";
    this->functionalGroupsIDsToWord[202] = "heterocyclic compound";
    this->functionalGroupsIDsToWord[203] = "alpha-aminoacid";
    this->functionalGroupsIDsToWord[204] = "alpha-hydroxyacid";

	this->fgsHierarchy.resize(205);
    this->fgsHierarchy[1] = -1;
    this->fgsHierarchy[2] = -1;
    this->fgsHierarchy[3] = -1;
    this->fgsHierarchy[4] = 3;
    this->fgsHierarchy[5] = 3;
    this->fgsHierarchy[6] = -1;
    this->fgsHierarchy[7] = 6;
    this->fgsHierarchy[8] = 6;
    this->fgsHierarchy[9] = -1;
    this->fgsHierarchy[10] = -1;
    this->fgsHierarchy[11] = -1;
    this->fgsHierarchy[12] = -1;
    this->fgsHierarchy[13] = -1;
    this->fgsHierarchy[14] = -1;
    this->fgsHierarchy[15] = -1;
    this->fgsHierarchy[16] = 15;
    this->fgsHierarchy[17] = -1;
    this->fgsHierarchy[18] = -1;
    this->fgsHierarchy[19] = -1;
    this->fgsHierarchy[20] = -1;
    this->fgsHierarchy[21] = -1;
    this->fgsHierarchy[22] = -1;
    this->fgsHierarchy[23] = -1;
    this->fgsHierarchy[24] = -1;
    this->fgsHierarchy[25] = -1;
    this->fgsHierarchy[26] = -1;
    this->fgsHierarchy[27] = -1;
    this->fgsHierarchy[28] = 27;
    this->fgsHierarchy[29] = 28;
    this->fgsHierarchy[30] = 28;
    this->fgsHierarchy[31] = 28;
    this->fgsHierarchy[32] = 30;
    this->fgsHierarchy[33] = 30;
    this->fgsHierarchy[34] = 27;
    this->fgsHierarchy[35] = 34;
    this->fgsHierarchy[36] = 27;
    this->fgsHierarchy[37] = -1;
    this->fgsHierarchy[38] = 37;
    this->fgsHierarchy[39] = 37;
    this->fgsHierarchy[40] = 37;
    this->fgsHierarchy[41] = -1;
    this->fgsHierarchy[42] = -1;
    this->fgsHierarchy[43] = -1;
    this->fgsHierarchy[44] = 43;
    this->fgsHierarchy[45] = -1;
    this->fgsHierarchy[46] = -1;
    this->fgsHierarchy[47] = -1;
    this->fgsHierarchy[48] = 47;
    this->fgsHierarchy[49] = 48;
    this->fgsHierarchy[50] = 48;
    this->fgsHierarchy[51] = 47;
    this->fgsHierarchy[52] = 51;
    this->fgsHierarchy[53] = 51;
    this->fgsHierarchy[54] = 51;
    this->fgsHierarchy[55] = 47;
    this->fgsHierarchy[56] = 55;
    this->fgsHierarchy[57] = 55;
    this->fgsHierarchy[58] = 55;
    this->fgsHierarchy[59] = -1;
    this->fgsHierarchy[60] = -1;
    this->fgsHierarchy[61] = -1;
    this->fgsHierarchy[62] = 61;
    this->fgsHierarchy[63] = 61;
    this->fgsHierarchy[64] = 61;
    this->fgsHierarchy[65] = 61;
    this->fgsHierarchy[66] = 61;
    this->fgsHierarchy[67] = 61;
    this->fgsHierarchy[68] = 61;
    this->fgsHierarchy[69] = 61;
    this->fgsHierarchy[70] = 61;
    this->fgsHierarchy[71] = 61;
    this->fgsHierarchy[72] = -1;
    this->fgsHierarchy[73] = 72;
    this->fgsHierarchy[74] = 72;
    this->fgsHierarchy[75] = -1;
    this->fgsHierarchy[76] = 75;
    this->fgsHierarchy[77] = 75;
    this->fgsHierarchy[78] = 75;
    this->fgsHierarchy[79] = 78;
    this->fgsHierarchy[80] = 75;
    this->fgsHierarchy[81] = 80;
    this->fgsHierarchy[82] = 80;
    this->fgsHierarchy[83] = 80;
	// TODO: check why 84 is not 80!
    this->fgsHierarchy[84] = -1;
    this->fgsHierarchy[85] = 75;
    this->fgsHierarchy[86] = 75;
    this->fgsHierarchy[87] = 75;
    this->fgsHierarchy[88] = 75;
    this->fgsHierarchy[89] = 75;
    this->fgsHierarchy[90] = -1;
    this->fgsHierarchy[91] = -1;
    this->fgsHierarchy[92] = 91;
    this->fgsHierarchy[93] = 91;
    this->fgsHierarchy[94] = 91;
    this->fgsHierarchy[95] = 91;
    this->fgsHierarchy[96] = -1;
    this->fgsHierarchy[97] = -1;
    this->fgsHierarchy[98] = -1;
    this->fgsHierarchy[99] = -1;
    this->fgsHierarchy[100] = 99;
    this->fgsHierarchy[101] = 99;
    this->fgsHierarchy[102] = 101;
    this->fgsHierarchy[103] = 99;
    this->fgsHierarchy[104] = 103;
    this->fgsHierarchy[105] = -1;
    this->fgsHierarchy[106] = -1;
    this->fgsHierarchy[107] = -1;
    this->fgsHierarchy[108] = -1;
    this->fgsHierarchy[109] = -1;
    this->fgsHierarchy[110] = 109;
    this->fgsHierarchy[111] = 109;
    this->fgsHierarchy[112] = -1;
    this->fgsHierarchy[113] = -1;
    this->fgsHierarchy[114] = 113;
    this->fgsHierarchy[115] = -1;
    this->fgsHierarchy[116] = -1;
    this->fgsHierarchy[117] = 116;
    this->fgsHierarchy[118] = 117;
    this->fgsHierarchy[119] = 117;
    this->fgsHierarchy[120] = 117;
    this->fgsHierarchy[121] = 117;
    this->fgsHierarchy[122] = 121;
    this->fgsHierarchy[123] = 121;
    this->fgsHierarchy[124] = 121;
    this->fgsHierarchy[125] = -1;
    this->fgsHierarchy[126] = 125;
    this->fgsHierarchy[127] = 125;
    this->fgsHierarchy[128] = 125;
    this->fgsHierarchy[129] = -1;
    this->fgsHierarchy[130] = 129;
    this->fgsHierarchy[131] = 129;
    this->fgsHierarchy[132] = 129;
    this->fgsHierarchy[133] = -1;
    this->fgsHierarchy[134] = -1;
    this->fgsHierarchy[135] = -1;
    this->fgsHierarchy[136] = -1;
    this->fgsHierarchy[137] = -1;
    this->fgsHierarchy[138] = -1;
    this->fgsHierarchy[139] = -1;
    this->fgsHierarchy[140] = -1;
    this->fgsHierarchy[141] = -1;
    this->fgsHierarchy[142] = -1;
    this->fgsHierarchy[143] = -1;
    this->fgsHierarchy[144] = -1;
    this->fgsHierarchy[145] = -1;
    this->fgsHierarchy[146] = -1;
    this->fgsHierarchy[147] = -1;
    this->fgsHierarchy[148] = -1;
    this->fgsHierarchy[149] = -1;
    this->fgsHierarchy[150] = -1;
    this->fgsHierarchy[151] = -1;
    this->fgsHierarchy[152] = -1;
    this->fgsHierarchy[153] = -1;
    this->fgsHierarchy[154] = 153;
    this->fgsHierarchy[155] = 153;
    this->fgsHierarchy[156] = 153;
    this->fgsHierarchy[157] = 153;
    this->fgsHierarchy[158] = 153;
    this->fgsHierarchy[159] = 153;
    this->fgsHierarchy[160] = 153;
    this->fgsHierarchy[161] = 153;
    this->fgsHierarchy[162] = 161;
    this->fgsHierarchy[163] = 161;
    this->fgsHierarchy[164] = -1;
    this->fgsHierarchy[165] = -1;
    this->fgsHierarchy[166] = -1;
    this->fgsHierarchy[167] = -1;
    this->fgsHierarchy[168] = -1;
    this->fgsHierarchy[169] = 168;
    this->fgsHierarchy[170] = 168;
    this->fgsHierarchy[171] = 168;
    this->fgsHierarchy[172] = 168;
    this->fgsHierarchy[173] = -1;
    this->fgsHierarchy[174] = 173;
    this->fgsHierarchy[175] = 173;
    this->fgsHierarchy[176] = 173;
    this->fgsHierarchy[177] = 173;
    this->fgsHierarchy[178] = -1;
    this->fgsHierarchy[179] = 178;
    this->fgsHierarchy[180] = 178;
    this->fgsHierarchy[181] = -1;
    this->fgsHierarchy[182] = 181;
    this->fgsHierarchy[183] = 181;
    this->fgsHierarchy[184] = 181;
    this->fgsHierarchy[185] = 181;
    this->fgsHierarchy[186] = -1;
    this->fgsHierarchy[187] = 186;
    this->fgsHierarchy[188] = 186;
    this->fgsHierarchy[189] = 186;
    this->fgsHierarchy[190] = 186;
    this->fgsHierarchy[191] = -1;
    this->fgsHierarchy[192] = 191;
    this->fgsHierarchy[193] = 191;
    this->fgsHierarchy[194] = -1;
    this->fgsHierarchy[195] = -1;
    this->fgsHierarchy[196] = -1;
    this->fgsHierarchy[197] = 196;
    this->fgsHierarchy[198] = 196;
    this->fgsHierarchy[199] = -1;
    this->fgsHierarchy[200] = -1;
    this->fgsHierarchy[201] = -1;
    this->fgsHierarchy[202] = -1;
    this->fgsHierarchy[203] = -1;
    this->fgsHierarchy[204] = -1;
   

}


/*
 * LigandModelCall::~LigandModelCall
 */
LigandModelCall::~LigandModelCall(void) {
}

bool LigandModelCall::SetCurrentLigandAndModel_byGlobalMdlID(unsigned int globalMdlID) { 
	if (globalMdlID <= (this->totalModelCount)) {
        this->currentGlobalMdlID = globalMdlID; 
		SetCurrentLigandAndModel(this->globalMdlID_to_ligID_mdlID[0][globalMdlID].GetX(),
            this->globalMdlID_to_ligID_mdlID[0][globalMdlID].GetY());
		return true;
	}
    return false;
}

bool LigandModelCall::SetCurrentLigandAndModel(unsigned int curLig, unsigned int curMdl) {
    if (this->ligandCount == 0) {
        this->currentLigand = 0;
        this->currentModel = 0;
        this->totalModelCount = 0;
        this->currentGlobalMdlID = -1;
        return false;
    }
    bool retval = true;
    this->currentLigand = curLig;
    this->currentModel = curMdl;
	// only for safety reasons (this conditions should never be true, but if true this is needed to avoid a crash)
    if (this->currentLigand >= this->ligandCount) {
        this->currentLigand = this->currentLigand % this->ligandCount;
        retval = false;
    }
    if (this->currentModel >= this->modelCount[this->currentLigand]) {
        this->currentModel = this->currentModel % this->modelCount[this->currentLigand];
        retval = false;
    }
    return retval;
}

void LigandModelCall::SetLigandAndModelCount(unsigned int ligCnt, const unsigned int* mdlCnt, unsigned int totalMdlCnt) { 
	this->ligandCount = ligCnt;
    this->modelCount = mdlCnt;
    this->totalModelCount = totalMdlCnt;
}

void LigandModelCall::SetTotalAtomPosCnt(unsigned int totalAtomPosCnt) {
	this->total_atomPosCnt = totalAtomPosCnt;
}


void LigandModelCall::SetBindigenergy(float bindingenergy) {
    this->currentBindingenrgy = bindingenergy;
}

void LigandModelCall::SetLigandname(vislib::StringA ligandname) {
    this->currentLigandname = ligandname;
}

void LigandModelCall::SetPDBQT_list_Filepath(vislib::StringA PDBQT_list_Filepath) {
    this->PDBQT_list_Filepath = PDBQT_list_Filepath;
}

void LigandModelCall::SetGlobalMdlID_to_ligID_mdlID(
    vislib::Array<vislib::math::Vector<uint, 2>>* globalMdlID_to_ligID_mdlID){
    this->globalMdlID_to_ligID_mdlID = globalMdlID_to_ligID_mdlID;
};

void LigandModelCall::SetSVGPath(std::string SVGPath) { 
	this->SVGPath = SVGPath;
}

void LigandModelCall::SetCurrentModelEfficiency(float modelEfficiency) {
    this->currentModelEfficiency = modelEfficiency;
}

void LigandModelCall::SetClusterAndCentroidData(
    vislib::Array<float>* modelCentroids,
	vislib::Array<float>* clusterCentroids,
	vislib::Array<int>* assignedClusters,
    vislib::Array<int>* clusterIndexSum, 
	vislib::Array<vislib::math::Cuboid<float>>* modelCentroidsCuboids,
    vislib::Array<vislib::math::Cuboid<float>>* clusterCentroidsCuboids,
    vislib::Array<vislib::math::Vector<int, 2>>* assign_modelID_sphereID, 
	vislib::Array<int> *clusterSizes,
    int minClusterSize,
	int maxClusterSize,
	int clusterDataHash) {
    this->c.modelCentroids = &modelCentroids[0];
    this->c.clusterCentroids = &clusterCentroids[0];
    this->c.assignedClusters = &assignedClusters[0];
    this->c.clusterIndexSum = &clusterIndexSum[0];
    this->c.modelCentroidsCuboids = &modelCentroidsCuboids[0];
    this->c.clusterCentroidsCuboids = &clusterCentroidsCuboids[0];
    this->c.assign_modelID_sphereID = &assign_modelID_sphereID[0];
    this->c.clusterSizes = &clusterSizes[0];
    this->c.minClusterSize = minClusterSize;
    this->c.maxClusterSize = maxClusterSize;
    this->c.maxClusterSize = maxClusterSize;
    this->c.clusterDataHash = clusterDataHash;
}

void LigandModelCall::SetInteractionForces_hbonds_protAcceptor(HProtAcceptors* hProtAcceptors) { 
	this->forces.hProtAcceptors = hProtAcceptors;
}

void LigandModelCall::SetInteractionForces_hbonds_protDonor(HProtDonors* hProtDonors) {
    this->forces.hProtDonors = hProtDonors;
}


void LigandModelCall::SetInteractionForces_hydrophobicInteractions(HydrophobicInteractions* hydrophobicInteractions){
    this->forces.hydrophobicInteractions = hydrophobicInteractions;
};

void LigandModelCall::SetInteractionForces_saltBridges(SaltBridges* saltBridges) {
    this->forces.saltBridges = saltBridges;
};
void LigandModelCall::SetInteractionForces_piStacks(PiStacks* piStacks) {
    this->forces.piStacks = piStacks;
};
void LigandModelCall::SetInteractionForces_piCationInteractions(PiCationInteractions* piCationInteractions) {
    this->forces.piCationInteractions = piCationInteractions;
};
void LigandModelCall::SetInteractionForces_halogenBonds(HalogenBonds* halogenBonds) {
    this->forces.halogenBonds = halogenBonds;
};
void LigandModelCall::SetInteractionForces_metalComplexes(MetalComplexes* metalComplexes) {
    this->forces.metalComplexes = metalComplexes;
};

void LigandModelCall::SetFunctionalGroupsOfLigand(FunctionalGroupsLigand* functionalGroupsLigand){
    this->functionalGroupsOfCurLigand = functionalGroupsLigand;
 };

void LigandModelCall::SetChemPropsOfLigand(ChemPropsLigand* chemPropsLigand){ 
	this->chemPropsOfCurLigand = chemPropsLigand;
 };

void LigandModelCall::SetDBSCANParams(float eps, int minPts, float minBindEnergy) { 
    this->c.DBSCANParams.eps = eps;
    this->c.DBSCANParams.minPts = minPts;
    this->c.DBSCANParams.minBindEnergy = minBindEnergy;
    this->c.DBSCANParams.paramsChanged = TRUE;
}

void LigandModelCall::SetlmcStatus(
    bool modelClustRenderer_busy, bool socketCommManager_busy, bool firstRenderCall_complete) {
    this->lmcStatus.isModelClustRendererBusy = modelClustRenderer_busy;
    this->lmcStatus.isSocketCommManagerBusy = socketCommManager_busy;
    this->lmcStatus.isFirstRenderCallComplete = firstRenderCall_complete;
}

void LigandModelCall::ResetDBSCANChanged() {
    this->c.DBSCANParams.paramsChanged = FALSE;
}

void LigandModelCall::SetWebData(WebData* webData) { 
	this->webData = webData; 
}

void LigandModelCall::SetClusteredFunctionalGroups(ClusteredFunctionalGroups* clusteredFunctionalGroups) {
    if (clusteredFunctionalGroups) {
        this->clusteredFunctionalGroups = clusteredFunctionalGroups;
        //} else if (this->clusteredFunctionalGroups->getDataHash() < clusteredFunctionalGroups->getDataHash()) {
    }
}

void LigandModelCall::ClusteredFunctionalGroups::setPointer_to_functionalGroupClusterSets(
    FGS_structN* functionalGroupClusterSets, int functionalGroupClusterSetsSize) {
    this->functionalGroupClusterSets = functionalGroupClusterSets;
    this->functionalGroupSetsCount = functionalGroupClusterSetsSize;
    this->dataHash++;
}




/* CLASS WEB DATA */
void LigandModelCall::WebData::UpdateWebData() {
	this->dataHash++; 
}
void LigandModelCall::WebData::UpdateGUIData() {
    this->GUIdataHash++;
    //UpdateWebData();
}
void LigandModelCall::WebData::UpdateSelectionData() {
    this->selection_dataHash++;
    //UpdateWebData();
}


