#include <limits>

enum instruction {div_HF};
enum CorL {Corl_C, Corl_L};

static float ComputeTable(instruction inst, CorL corl){
    switch (inst)
    {
    case div_HF:
        if(corl == Corl_C){
            return 15;
        }else{
            return 5;
        }
        break;
    default:
        return std::numeric_limits<float>::infinity();
        break;
    }
}