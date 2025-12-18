#include "convergenceNode.h"

#include "autoTarget.h"
#include "targetError.h"

#include "readInput.h"
#include "printOutput.h"
#include "farm.h"
#include "asciiFile.h"
#include "metric.h"

#include "autoconverge.h"

const std::shared_ptr<AutoTarget> ConvergenceNode::target() const {return _root->target();}
std::shared_ptr<Metric> ConvergenceNode::metric() const {return _root->metric();}
std::string ConvergenceNode::reference() const {return _root->reference();}
std::string ConvergenceNode::referenceRoot() const  {return _root->reference();}


const std::string ConvergenceNode::macro="_MACRO_";
std::string ConvergenceNode::estimate="EstimateError";

static double goodRatio=1e-1; // ratio between subsequent differences that is considered good for fitting
static double minRatio=goodRatio*2; // acceptable ratio for fitting
static double distMax=0.3; // avoid using with more than this error
static int maxStepsUp=9;   // default maximal steps upward from initial
static int maxStepsDown=2; // default maximal steps downward from initial

static std::string name(std::string Entry){
    size_t p=Entry.find(':');
    return p==std::string::npos?Entry:tools::cropString(Entry.substr(p+1));
}

std::string ConvergenceNode::strNode(int Level) const{
    if(Level!=Tree_defaultKind)ABORT("only default");
    Str s(""," ");
    s=s+_nam+" "+_cat+_parMin+_parMax;
    return std::move(s);
}

ConvergenceNode::ConvergenceNode(ReadInput &Inp, const Autoconverge* Root, int line):_size(0),_root(Root)
{
    std::string preCat(macro),atLine,strInit,strStep;
    Inp.read("Autoconverge","category",_cat,preCat,"Category of input item to vary (blank - use previous) or \""+macro+"\" for varying macros",line);
    if(_cat!=preCat)preCat=_cat;
    if(tools::toLower(_cat).find("macro")<4 and _cat!=macro)ABORT("for varying macros use Autoconverge: category  \""+macro+"_\", found: "+_cat);
    Inp.read("Autoconverge","name",_nam,"","name of input item to vary",line);
    Inp.read("Autoconverge","atLine",atLine,ReadInput::anyName,"enter matchName=val to vary entry at line where matchName==val in category",line);
    Inp.read("Autoconverge","step",strStep,"0.2 of initial","initial step for varying parameter, specify as \"10\" or \"0.3 of initial\"",line);
    Inp.read("Autoconverge","lowerLimit",_parMin,tools::str( DBL_MAX),Str("","")+"lower limit for parameter (default: present-"+maxStepsDown+"*step)",line);
    Inp.read("Autoconverge","upperLimit",_parMax,tools::str(-DBL_MAX),Str("","")+"upper limit for parameter (default: present+"+maxStepsUp+"*step)",line);

    if(Inp.found("Autoconverge")){
        // construct data on node
        if(_cat!=macro and not Inp.found(_cat,_nam))ABORT("Autoconverge entry "+_cat+": "+_nam+" not found in input");

        // get input as a matrix, return in _row,_col location of present parameter in matrix
        _inputRows=inputRowCol(Inp,atLine,_row,_col);

        if(_col<0 or _row<0)
            ABORT("Autoconverge: "+_cat+ReadInput::categoryTerminator+" "+_nam+" atLine "+atLine+" could not be located input file");
        _pars.assign(1,tools::string_to_double(inputRows()[_row][_col]));


        if(strStep.find("of initial")!=std::string::npos){
            _parStep=std::abs(_pars.front())*tools::string_to_double(strStep.substr(0,strStep.find("of initial")+1));
        }
        else
            _parStep=tools::string_to_double(strStep);
        if(_parStep<0)ABORT("must have postive increment")
                if(_parStep>1)_parStep=int(_parStep);

        if(_parMax<-DBL_MAX/2)_parMax=_pars.front()+maxStepsUp*_parStep;
        if(_parMin> DBL_MAX/2)_parMin=_pars.front()-maxStepsDown*_parStep;
        _size=2; // will compute 2 additional results
    }

}

int ConvergenceNode::parameterRow() const{
    for(size_t k=0;k<_inputRows.size();k++){
        std::vector<std::string> Row=_inputRows[k];
        if(Row.size()){
            if(_cat==macro and Row[0]=="#define" and Row[1].find(_nam)!=std::string::npos)return k;
            if(Row[0].find(_cat+ReadInput::categoryTerminator)==0
                    and Row[0].substr(Row[0].find(ReadInput::categoryTerminator)).find(_nam)!=std::string::npos)return k;
        }
    }
    return _row+1; // parameterRow after value indicates not found
}


void ConvergenceNode::getDeltas(std::vector<double>& CPars, std::vector<double> &Delta) const
{
    std::vector<double> IPars(_pars);
    std::vector<std::string> Inps(_inps);
    IPars.resize(Inps.size());
    tools::sortByKey(IPars,Inps);

    // determine convergence factors between subsequent and to largest parameter

    // check for completed runs:
    std::vector<size_t> comp;
    CPars.clear();
    for(size_t k=0;k<Inps.size();k++){
        if(target()->cols(Inps[k]).size()>0){
            comp.push_back(k);
            CPars.push_back(IPars[k]);
        }
    }

    Delta.clear();
    for(size_t k=1;k<comp.size();k++)
        Delta.push_back(distance(Inps[comp.back()],Inps[comp[k-1]]));

}

void ConvergenceNode::parametersAndErrors(std::vector<double> &Pars, std::vector<double> &Estims) const{
    getDeltas(Pars,Estims);
    Estims=errorFunction()->estimates(Pars,Estims);
}

std::string ConvergenceNode::writeConvergence(){
    std::vector<double>cpars,delta;
    getDeltas(cpars,delta);

    AsciiFile f(tools::stringStrip(reference(),"inpc")+"convData");
    std::vector<std::string> head;
    std::string inps;
    for(size_t k=0;k<cpars.size();k++)inps+="  "+_inps[k];
    head.push_back(inps);
    head.push_back(_nam+"         d_best      d_prev    d_prev/d_prev-1");
    f.writeComments(head);
    f.writeCols({{cpars.begin(),cpars.end()-1},delta});

    return f.name();
}

double ConvergenceNode::distance(std::string DirA, std::string DirB) const {

    std::vector<std::vector<double>> colsA(target()->cols(DirA));
    std::vector<std::vector<double>> colsB(target()->cols(DirB));

    double dsqu=0;
    for(size_t c=0;c<colsA.size();c++){
        dsqu+=std::pow(metric()->distance(colsA[c],colsB[c]),2);
    }
    return sqrt(dsqu/colsA[0].size());
}

// Algorithm
// - find error difference that best approximates desired error
// - locate two differences in the vicinity that show improvement by a "goodRatio" at least
// - if none, take best approximation
// - limit parameter increase by maxStep
double ConvergenceNode::nextParameter(std::vector<double> &CPars, std::vector<double> &EstimErrs, double &ParAtTarget) const {
    /// less than 3 calculations started - no meaningful error estimate possible
    if(_pars.size()<3)DEVABORT("not supposed to call this")

            // wait until at least 3 calcs are completed
            std::vector<double> Deltas,alphas,gammas;
    do {
        getDeltas(CPars,Deltas);
        if(CPars.size()<_inps.size())tools::sleepSeconds(0.1);
    } while (CPars.size()<3 or CPars.size()+3<_inps.size());

    if(not std::is_sorted(CPars.begin(),CPars.end()))DEVABORT(Sstr+"need parameters sorted by increasing size, got: "+CPars);

    if(Deltas.size()==0)return _parMin;

    //-------- determine close error and vicinity [klo,kup] with sufficient variation ---------
    size_t ref=CPars.size(),kup,klo;

    // find reference among sorted parameters
    ref=std::find(CPars.begin(),CPars.end(),_pars[0])-CPars.begin();

    // ascend from ref to larger until good ratio (or end)
    for(kup=ref;kup+1<Deltas.size() and Deltas[kup]/Deltas[ref]>goodRatio;kup++);

    // descend from ref to lower parameter until good ratio (or beginning)
    for(klo=ref;klo>0 and Deltas[ref]/Deltas[klo]>goodRatio;klo--);

    if(klo==kup)kup=klo+1;

    // must not use last available par for error ratio
    if(kup+1==CPars.size()){
        kup--;
        if(klo>0)klo--;
    }

    // avoid estimates from large errors (above distMax)
    while(klo<Deltas.size()-1 and Deltas[klo]>distMax){
        klo=std::min(klo+1,Deltas.size()-1);
        kup=std::min(kup+1,Deltas.size()-1);
    }

    double next;
    if(klo<kup){
        //---- fit TargetError to Deltas -----------------------------------------------------
        errorFunction()=TargetError::factory(ref,klo,kup,_parMin,_parMax,CPars,Deltas,0.1);
        EstimErrs.clear();
        for(double p: CPars)EstimErrs.push_back(errorFunction()->err(p));
        // do not extrapolate to large values, use Deltas as error estimates instead
        for(size_t k=0;k<Deltas.size();k++)if(Deltas[k]>EstimErrs.back()*10)EstimErrs[k]=Deltas[k];
    }
    else {
        errorFunction().reset(new ErrorNoise(CPars.back(),_parMax));
        next=CPars.back()+std::max(1.,CPars.back()*0.2);
        EstimErrs=Deltas;
        EstimErrs.push_back(0);
    }

    // where differences to best are much larger than last increment
    // use increment for error estimate, not extrapolation of fit
    for(size_t k=0;k<Deltas.size();k++)
        if(Deltas.back()*100<Deltas[k])EstimErrs[k]=Deltas[k];

    // suggest next parameter
    next=_pars[0]; // we have 3, try to accept


    if(Deltas[kup]/Deltas[klo]>minRatio){
        // not a sane fit - try to span two "good ratios" from closest down to smaller parameters

        next=std::max(_parMin,std::min(_parMax,errorFunction()->inv(*std::min_element(EstimErrs.begin(),EstimErrs.end())/std::pow(goodRatio,1))));
        next=std::min(next,_parMax);

        if(next<_parMin){
            // below lower limit, extend upward instead
            next=std::max(_parMin,std::min(_parMax,errorFunction()->inv(errorFunction()->err(_parMin)*std::pow(goodRatio,2))));
            // new span already tried, increase overall
            while(next<CPars.back() and next<_parMax)next=CPars.back()+_parStep;
        }

        // avoid too close repetition
        while(std::count_if(_pars.begin(),_pars.end(),
                            [&](double val){return std::abs(val-next)<_parStep; } )
              and next<_parMax){
            next+=_parStep;
        }
    }

    // avoid too large single steps
    double maxStep=std::min(_parStep*3,0.3*(_parMax-_parMin));
    next=std::min(int(next),int(CPars.back()+maxStep));

    ParAtTarget=_pars[0];

    return next;

}


std::vector<std::vector<std::string>> ConvergenceNode::inputRowCol(ReadInput& Inp, std::string AtLine, int &Row, int &Col){
    std::vector<std::vector<std::string>> rowsCols;
    if(AtLine!=ReadInput::anyName){
        std::vector<std::string> inpLines=Inp.inputLines();
        std::string linNam=AtLine.substr(0,AtLine.find_first_of("=["));
        bool inCategory=false;
        int rowAfterCat=-1;
        for(auto l: inpLines){
            if(l.find("#define ")==0){
                rowsCols.push_back(tools::splitString(l,' '));
                if(_cat==macro and _nam==rowsCols.back()[1]){
                    Col=2;
                    Row=rowsCols.size()-1;
                    inCategory=true;
                }
            }
            else
                rowsCols.push_back(tools::splitString(l,',',ReadInput::leftBracket+ReadInput::quote,ReadInput::rightBracket+ReadInput::quote));
            for(auto &r: rowsCols.back())r=tools::cropString(r);
            if(ReadInput::isCategoryLine(l) and l.find(_cat+ReadInput::categoryTerminator)==0){
                Col=0;
                while(name(rowsCols.back()[Col])!=_nam)Col++;
                inCategory=true;
            }
            if(inCategory){
                rowAfterCat++;
                // determine row
                if(AtLine==linNam){
                    if(rowAfterCat==tools::string_to_int(linNam))
                        Row=rowsCols.size()-1;
                }
                else {
                    if(_nam==linNam and tools::string_to_int(tools::stringInBetween(AtLine,"[","]")))
                        Row=rowsCols.size()-1;
                }
                if(Inp.endCategory(_cat,rowAfterCat))inCategory=false;
            }
        }
    }
    return rowsCols;
}

















