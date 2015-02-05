ModelMatrixDimension := 64;

// -------------------- initialize two matrix -------------------------------

function PopulateModelMatrix (ModelMatrixName&, pi, mu) {
	if (!ModelMatrixDimension)
        {
             ModelMatrixDimension = 64;
        }
       
    ModelMatrixName = {ModelMatrixDimension,ModelMatrixDimension};
	
	modelDefString = "";
	// 4*demension*dimension
    modelDefString*65536;
	for (h=0; h<ModelMatrixDimension; h=h+1) {
		for (v = 0; v<ModelMatrixDimension; v=v+1) {
			if(h < 100 && v < 100) {
				if(h == v+1 || h == v-1) {
					modelDefString*("ModelMatrixName["+h+"]["+v+"] := " + "(pi" + "+" + "mu)*t1" + ";\n");
				} else {
					if(h == v){
						if(h == 0 || h == ModelMatrixDimension-1) {
							// -(pi*(ModelMatrixDimension-1) + mu)
							modelDefString*("ModelMatrixName["+h+"]["+v+"] := 0-(" + "pi" + "*(" + ModelMatrixDimension + "-" + 1 + ")+" + "mu" + ")*t1;\n");
						} else {
							// -(pi*(ModelMatrixDimension-1) + 2*mu)
							modelDefString*("ModelMatrixName["+h+"]["+v+"] := 0-(" + "pi" + "*(" + ModelMatrixDimension + "-" + 1 + ")+2*" + "mu" + ")*t1;\n");
						}
					} else{
						modelDefString*("ModelMatrixName["+h+"]["+v+"] :=" + "pi*t1" + ";\n");
					}
				}
			} else {
				modelDefString*("ModelMatrixName["+h+"]["+v+"] :=" + "0" + ";\n");
			}
		}
	}
	modelDefString*0;
    ExecuteCommands (modelDefString);

    return 0;
}

function PopulateModelMatrixFreq (ModelMatrixName&) {
	if (!ModelMatrixDimension)
        {
             ModelMatrixDimension = 64;
        }
       
    ModelMatrixName = {1,ModelMatrixDimension};
	
	modelDefString = "";
    modelDefString*512;
	for (h=0; h<ModelMatrixDimension; h=h+1) {
		if(h < 100) {
			modelDefString*("ModelMatrixName["+0+"]["+h+"] := " + 1/ModelMatrixDimension+ ";\n");
		} else {
			modelDefString*("ModelMatrixName["+0+"]["+h+"] := " + 0 + ";\n");
		}
	}
	modelDefString*0;
    ExecuteCommands (modelDefString);
    return 0;
} 


// ----------------------constructing the matrix of logistic regression with only first two columns have values ------------------------

function ConstructLGProbMatrix(ModelMatrixName&, a, b) {
	ModelMatrixName = {ModelMatrixDimension,ModelMatrixDimension};
	
	modelDefString = "";
	// 4*demension*dimension
    modelDefString*65536;
	for(h=0; h<ModelMatrixDimension; h=h+1) {
		if(h < 100) {
			modelDefString*("ModelMatrixName["+h+"]["+0+"] := " + "1 - 1/(1+Exp(-" + "a" + "-" + "b"+ "*" + h/ModelMatrixDimension + "))"  + ";\n");  // P(Y=0)
			modelDefString*("ModelMatrixName["+h+"]["+1+"] := " + "1/(1+Exp(-" + "a" + "-" + "b"+ "*" + h/ModelMatrixDimension + "))"  + ";\n");  // P(Y=1)
			for(s=2; s<ModelMatrixDimension; s=s+1) {
				modelDefString*("ModelMatrixName["+h+"]["+s+"] := " + 0 + ";\n");
			}
		} else {
			for(s=0;s<ModelMatrixDimension; s=s+1) {
				modelDefString*("ModelMatrixName["+h+"]["+s+"] := " + 0 + ";\n");
			}
		}
	}
	modelDefString*0;
    ExecuteCommands (modelDefString);
    return 0;

}


// ---------------Sergei's code -------------------------

function define_logistic_mapper (dimension, model_name, logistic_regression_parameter_prefix) {
    // see http://en.wikipedia.org/wiki/Logistic_regression
    
    // zero NxN matrix   
    define_logistic_mapper.matrix = model_name + ".transition_matrix";    
    ExecuteCommands ("`define_logistic_mapper.matrix` = {dimension, dimension};");

    // fake Nx1 frequencies matrix [needed to define the model, but will be ignored]   
    define_logistic_mapper.frequencies = model_name + ".frequencies";
    ExecuteCommands ("`define_logistic_mapper.frequencies` = {dimension,1}[\"1/dimension\"]"); // set all entries to 1/dimension
    
    // logistic model parameters    
    define_logistic_mapper.intercept = logistic_regression_parameter_prefix + ".beta0";
    define_logistic_mapper.slope     = logistic_regression_parameter_prefix + ".beta1";
    
    ExecuteCommands ("global `define_logistic_mapper.intercept` = 0.0");
    *define_logistic_mapper.intercept :> -1e10;
    *define_logistic_mapper.intercept :< 1e10;

    ExecuteCommands ("global `define_logistic_mapper.slope` = 0.1");
    *define_logistic_mapper.slope :> -1e10;
    *define_logistic_mapper.slope :< 1e10;
    
    for ( define_logistic_mapper.index = 0;  define_logistic_mapper.index < dimension; define_logistic_mapper.index += 1) {
    
        ExecuteCommands ("`define_logistic_mapper.matrix`[define_logistic_mapper.index][define_logistic_mapper.index] = 0");
        
        ExecuteCommands ("`define_logistic_mapper.matrix`[define_logistic_mapper.index][1] := " + 
                "1/(1+Exp (-`define_logistic_mapper.intercept` - `define_logistic_mapper.slope` * (" + (define_logistic_mapper.index - dimension$2) + ")))") ;
        
        // probability of positive phenotype given "hidden value = index"
 
        ExecuteCommands ("`define_logistic_mapper.matrix`[define_logistic_mapper.index][0] := " + 
                "1- 1/(1+Exp (-`define_logistic_mapper.intercept` - `define_logistic_mapper.slope` * (" + (define_logistic_mapper.index - dimension$2) + ")))") ;
            
        // 1 - probability of positive phenotype given "hidden value = index"
    }
    fprintf(stdout,"Model `model_name` = (\"`define_logistic_mapper.matrix`\", `define_logistic_mapper.frequencies`, EXPLICIT_FORM_MATRIX_EXPONENTIAL)");
    ExecuteCommands ("Model `model_name` = (\"`define_logistic_mapper.matrix`\", `define_logistic_mapper.frequencies`, EXPLICIT_FORM_MATRIX_EXPONENTIAL)");
    
}

function tree_extender (tree_id, leaf_suffix, model_name, original_model_name) {
    tree_extender.leaf_count = TipCount (*tree_id);
    for (tree_extender.counter = 0; tree_extender.counter < tree_extender.leaf_count; tree_extender.counter += 1) {
        tree_extender.leaf_name = TipName (*tree_id, tree_extender.counter);
        
        //fprintf(stdout, (*tree_id));
		//fprintf(stdout, {"WHERE": tree_extender.leaf_name, "NAME": tree_extender.leaf_name + leaf_suffix});
		//fprintf(stdout, "\n");
		// ExecuteCommands ("ReplicateConstraint(this1." + tree_extender.leaf_name + ".t1:=this2." + tree_extender.leaf_name+ ".t, newTree, givenTree)");
        (*tree_id) + {"WHERE": tree_extender.leaf_name,
                      "NAME": tree_extender.leaf_name + leaf_suffix};
        //fprintf(stdout,"SetParameter (" + tree_id + "." + tree_extender.leaf_name + leaf_suffix + ", MODEL, `model_name`)");
        ExecuteCommands ("SetParameter (" + tree_id + "." + tree_extender.leaf_name + leaf_suffix + ", MODEL, `model_name`)");
        //ExecuteCommands ("SetParameter (" + tree_id + "." + tree_extender.leaf_name + ", MODEL, `original_model_name`)");
		//ExecuteCommands ("ReplicateConstraint(this1." + tree_extender.leaf_name + ".t1:=this2." + tree_extender.leaf_name+ ".t, newTree, givenTree)");
    }
}

// --------------- end of Sergei's code----------------



// starting the main code

// input the observations and treeString

SetDialogPrompt ("Specify the data file"); 
fprintf(stdout, "Please input the data file\n");
DataSet myData = ReadDataFile (PROMPT_FOR_FILE);
fprintf (stdout, "Loaded ", myData.species, " sequences with ", myData.sites, " sites from ",LAST_FILE_PATH,"\n");

SetDialogPrompt ("Specify the phenotype file"); 
fprintf(stdout, "Please input the phenotype file\n");
DataSet myData2 = ReadDataFile (PROMPT_FOR_FILE);
fprintf (stdout, "Loaded observations of ", myData2.species, " sequences", " from ",LAST_FILE_PATH,"\n");

//DataSetFilter myFilter1 = CreateFilter (myData, 1);
//DataSetFilter myFilter2 = CreateFilter (myData2, 2);

SetDialogPrompt ("Specify the tree file"); 
fprintf(stdout, "Please input the tree string\n");
fscanf (PROMPT_FOR_FILE, "Raw", treeString);

fprintf(stdout, "The tree is like the following :\n" + treeString + "\n");

fprintf(stdout, "The input part finishes\n");


/*--Code--*/
//

/*---------Estimate Branch Lengths Using Nucleotide Model-----------------------*/
DataSetFilter myFilter = CreateFilter (myData,1);
HarvestFrequencies (obsNucFreqs, myFilter, 1, 1, 1);


/*---------Begin Setting up custom nuc model---------------*/
/*The parameters that the nucModelString indexes into*/
nucBiasMult = {{"AC*","","AT*","CG*","CT*","GT*"}};
nucModelString 	= "012345";

/*Initialises and populates a matrix of string values nuc multipliers*/
customRateString = {{"*","","",""}
					{"","*","",""}
					{"","","*",""}
					{"","","","*"}
				 };
				 
for (i=0; i<3; i+=1)
{
	shift = -(i == 0);
	for(j=i+1;j<4;j += 1)
	{
		customRateString[i][j] = Eval ("nucBiasMult["+nucModelString[j+i+shift]+"]");
		customRateString[j][i] = customRateString[i][j];
	}
}

global AC = 1; global AT = 1; global CG = 1; global CT = 1; global GT = 1; 

/*To set up a nucleotide model, populate a string version of the rate matrix*/
modelDefString = "";
modelDefString * 16384;

modelDefString* "{" ;
for (i=0; i<4; i += 1)
{
	modelDefString*("{");
	for(j=0;j<4;j=j+1)
	{
		if(j>0)
		{
			modelDefString*(",");
		}
		if(i==j)
		{
			modelDefString*("*");
		}
		else
		{
			modelDefString*(customRateString[i][j]+"t");
		}
	}
	modelDefString*"}";
}
modelDefString*"}";
modelDefString*0;
ExecuteCommands("nucModel = " + modelDefString);
Model RT = (nucModel, obsNucFreqs);

//fprintf(stdout, "before the nuc model fit the nucModel matrix:" + nucModel);


/*-----------End of defining nuc model------------*/


/*---------start nuc model fit----------*/
Tree tempTree = treeString;
treeString = RerootTree(treeString,BranchName(tempTree,0));

Tree givenTree = treeString;

fprintf (stdout, "\n\n[PHASE 1. Estimating Branch Lengths using a Nucleotide Model]\n");
LikelihoodFunction theLikFun = (myFilter, givenTree, obsFreqs);
Optimize (paramValues, theLikFun);
fprintf (stdout, theLikFun);
fprintf (stdout, "\n", "The end of nuc model");
fprintf (stdout, "\n");
fprintf (stdout, "\n\n[PHASE 1 finished]\n");
/*---------end nuc model fit----------*/

/*--------Forever constrain nuc rates-------*/
AC := AC__; AT := AT__; CG := CG__; CT := CT__; GT := GT__;


/*--------- creating transition matrix --------------*/
// -------------------- define params -------------
// parameters for Q matrix
global mu  = 1;
mu :< 10;
global pi = 1;
pi := 0;
// parameters for logistic regression
global a = -0.5;
ClearConstraints(a);
a:>-10; a:<10;
global b = ModelMatrixDimension;
b:< 100;

//fprintf(stdout, "\n");


// --------------- initialize matrix -------------------
PopulateModelMatrix("Q", pi, mu);
fprintf(stdout, Q);
PopulateModelMatrixFreq("eqFrequencies");
fprintf(stdout, eqFrequencies);
ConstructLGProbMatrix("LGProb",a,b);
fprintf(stdout,LGProb);
fprintf(stdout, "\n");
/*--------- end of creating transition matrix ----------*/



/* --------- extend the giventree by tree_tender function ---------*/
fprintf (stdout, "\n\n[PHASE 2. Estimating parameters for transition matrix and logistic regression parameters]\n");
Model m1 = (Q, eqFrequencies);
Model LRModel = ("LGProb", eqFrequencies, EXPLICIT_FORM_MATRIX_EXPONENTIAL);
INCLUDE_MODEL_SPECS = 1;
UseModel (USE_NO_MODEL);
newTreeString = "" + givenTree;
newTreeString = newTreeString^{{"{RT}"}{"{m1}"}};
Tree newTree = "" + newTreeString;


/* set constraints */
// ClearConstraints    (newTree);
ReplicateConstraint ("this1.?.t1:=this2.?.t__",newTree,givenTree);
//fprintf (stdout, givenTree.CY002120_1998_12_New_York_254.t, "\n", newTree.CY002120_1998_12_New_York_254.mu, "\n", newTree.CY002120_1998_12_New_York_254.t1);


fprintf (stdout, "\nBEFORE\n", newTree, "\n");
tree_extender ("newTree", "_phenotype", "LRModel", "m1");
fprintf (stdout, "\nAFTER\n", newTree, "\n");
/*------- end of tree_extender -------------*/

//fprintf (stdout, givenTree.CY002120_1998_12_New_York_254.t, " ", newTree.CY002120_1998_12_New_York_254.mu, " ", newTree.CY002120_1998_12_New_York_254.t1, "\n");
//fprintf (stdout, newTree.CY002120_1998_12_New_York_254.t1, "  ", newTree.CY002120_1998_12_New_York_254.mu, " ", newTree.CY002120_1998_12_New_York_254.pi,"\n");
//fprintf (stdout, newTree.CY002120_1998_12_New_York_254_phenotype.t1, "  ", newTree.CY002120_1998_12_New_York_254_phenotype.a, " ", newTree.CY002120_1998_12_New_York_254_phenotype.b,"\n");
//fprintf (stdout, newTree.Node32.t1, "  ", newTree.Node32.mu, " ", newTree.Node32.pi,"\n");


/* --------- estimate paramaters for transition matrix and Logistic regression matrix -------*/

UseModel(m1);
DataSetFilter myFilter2 = CreateFilter (myData2,6);
LikelihoodFunction lf2 = (myFilter2, newTree);
Optimize (myRes, lf2);
fprintf (stdout, lf2);

summaryPath := "K:\\winter2015\\project\\simple_test_data\\summary"; 
Export (export_string, lf2);
fitPath=summaryPath+".fit";
fprintf(fitPath,export_string);


//fprintf (stdout, myRes);
fprintf (stdout, "\n\n[PHASE 2 finished]\n");

//fprintf (stdout, givenTree.CY002120_1998_12_New_York_254.t, " ", givenTree.CY002120_1998_12_New_York_254_phenotype.mu, " ", newTree.CY002120_1998_12_New_York_254.t1, "\n");
//fprintf (stdout, newTree.CY002120_1998_12_New_York_254.t1, "  ", newTree.CY002120_1998_12_New_York_254.mu, " ", newTree.CY002120_1998_12_New_York_254.pi,"\n");
//fprintf (stdout, newTree.CY002120_1998_12_New_York_254_phenotype.t1, "  ", newTree.CY002120_1998_12_New_York_254_phenotype.a, " ", newTree.CY002120_1998_12_New_York_254_phenotype.b,"\n");
//fprintf (stdout, newTree.Node32.t1, "  ", newTree.Node32.mu, " ", newTree.Node32.pi,"\n");
