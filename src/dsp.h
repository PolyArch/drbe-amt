#include <limits>

enum instruction {Div_HF, Div, Sqrt, Arctan, Sin, Cos, Cosi, division, mod};
enum CorL {Corl_C, Corl_L};

static float ComputeTable(instruction inst, CorL corl){
    switch (inst)
    {
    case Div_HF:
        if(corl == Corl_C){
            return 15;
        }else{
            return 5;
        }
        break;
    case Div:
        if(corl == Corl_C){
            return 10;
        }else{
            return 4;
        }
        break;
    case Sqrt:
        if(corl == Corl_C){
            return 10;
        }else{
            return 5;
        }
        break;
    case Arctan:
        if(corl == Corl_C){
            return 24;
        }else{
            return 5;
        }
        break;
    case Sin:
        if(corl == Corl_C){
            return 10;
        }else{
            return 5;
        }
        break;
    case Cos:
        if(corl == Corl_C){
            return 10;
        }else{
            return 5;
        }
        break;
    case Cosi:
        if(corl == Corl_C){
            return 20;
        }else{
            return 6;
        }
        break;
    case division:
        return 20;
    case mod:
        return 1;
    default:
        return std::numeric_limits<float>::infinity();
        break;
    }
}