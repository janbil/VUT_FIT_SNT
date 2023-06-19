#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <math.h>
#include <experimental/filesystem>
#define E 2.71828182846
using std::experimental::filesystem::directory_iterator;
//Struktura reprezentujici aktivitu
struct activity {
    int id;
    int duration;         
    int resource_needed [4];
    int num_succ;
    int succs[31];
}; 
//Globalni promenne pro maximalni dostupny pocet zdroju a pocet aktivit v projektu
int resource_avaiabilities[4];
int num_of_act;

//Trida, ktera vypocita sousedni reseni a jeho makespan
class SA {
  public:
    int resource_avaiability[4];
    activity activities [62]; //Potreba zmenit pro zmenu datasetu J30 = [32], J60 = [62]
    std::vector<int> given_list_repr;
    std::vector<int> list_repr;
    std::vector<int> times;
    std::vector<int> end_times;
    std::set<int> currently_active;
    std::set<int> available_acts;
    int curr_index;
    int curr_time;
    int curr_available_res[4];
    //generuje sousedni reseni
    void generate_neighbor(){
        int rand_act = rand() % (num_of_act-2) + 1;
        int latest_pre = latest_pred(rand_act);
        int lower_bound;
        for(int i =0; i<given_list_repr.size(); i++){
            if(given_list_repr[i] == latest_pre){
                lower_bound = i;
            }
        }
        int rand_pos;
        for(int i =0; i<given_list_repr.size(); i++){
            if(given_list_repr[i] == rand_act){
                rand_pos = i;
            }
        }
        int earliest = earliest_successor(rand_act);
        int new_pos = rand() % (earliest-lower_bound-1) + lower_bound+1;
        int iter = 0;
        while(new_pos == rand_pos){
            new_pos = rand() % (earliest-lower_bound-1) + lower_bound+1;
            if(iter>5){
                generate_neighbor();
                return;
            }
            iter++;
        }
        if(new_pos > rand_pos){
            int to_new_pos = given_list_repr[rand_pos];
            for(int i = rand_pos+1; i<= new_pos; i++){
                given_list_repr[i-1] = given_list_repr[i];
            }
            given_list_repr[new_pos] = to_new_pos;
        } else {
            int to_new_pos = given_list_repr[rand_pos];
            for(int i = rand_pos; i> new_pos; i--){
                given_list_repr[i] = given_list_repr[i-1];
            }
            given_list_repr[new_pos] = to_new_pos;
        }
    }
    //Kontroluje, zda se nevygenerovalo reseni porusujici souslednost. V podstate debugging funkce
    bool check_constraints(std::vector<int> new_list){
        int curr_ind = 0;
        for(auto act : new_list){
            for(auto pre : pred(activities[act])){
                bool found = false;
                for(int i=0; i<curr_ind; i++){
                    if(new_list[i] == pre){
                        found = true;
                    }
                }
                if(!found){
                    std::cout<<"Chyba pred "<<pre << " act: " << act <<std::endl;
                    return false;
                }
            }
            for(int ii = 0; ii<activities[act].num_succ; ii++){
                bool found = false;
                for(int j=curr_ind; j<new_list.size(); j++){
                    if(new_list[j] == activities[act].succs[ii]){
                        found = true;
                    }
                }
                if(!found){
                    std::cout<<"Chyba succ "<<std::endl;
                    return false;
                }
            }
            curr_ind++;
        }
        return true;
    }

    //Najde nejblizsiho naslednika
    int earliest_successor(int act){
        int earliest = INT32_MAX;
        for(int i =0; i< given_list_repr.size() ; i++){
            for(int j =0;j<activities[act].num_succ;j++){
                if(given_list_repr[i] == activities[act].succs[j]){
                    if(i < earliest){
                        
                        earliest = i;
                    }
                }
            }
        }
        return earliest;
    }

    //Najde nejposlednejsiho predchudce
    int latest_pred(int act){
        int latest_index=0;
        for(int i = 0; i<given_list_repr.size() ; i++){
            for(auto pred : pred(activities[act])){
                if(given_list_repr[i] == pred){
                    if(i>latest_index){
                        latest_index = i;
                    }
                }
            }
        }
        return given_list_repr[latest_index];
    }

    //Hlavni funkce, nastavi parametry, najde vedlejsi reseni a odsimuluje aby nasel makespan.
    void compute_timespan (){
        init();
        generate_neighbor();
        if(!check_constraints(given_list_repr)){
            std::cout<<"chybka"<<std::endl;
            exit(1);
        }
        while(list_repr.size()<num_of_act){
            step();
        }
    }

  private:  
    //Inicializace parametru a naplanovani prvni dummy aktivity.
    void init (){
        
        curr_index = 0;
        curr_time=0;
        for(int i=0;i<4;i++){
            curr_available_res[i] = resource_avaiability[i];
        }
        available_acts.insert(0);
        schedule_activity(activities[0]);
    }

    //Krok simulace
    void step(){
        int act_to_plan = -1;
        
        for(auto ac : available_acts){
            if(ac == given_list_repr[curr_index]){
                act_to_plan = given_list_repr[curr_index];
            }
        }
        //Nemuzu naplanovat spravnou aktivitu, protoze bych porusil constraints, 
        //musim dokoncit nejblizsi rozdelanou a posunout cas
        if(act_to_plan == -1){
            if(currently_active.size() == 0){
                std::cout<<"Error" << std::endl;
            }
            
            int lowest_time = INT32_MAX;
            for(auto time : end_times){
                if(time > curr_time){
                    if(time<lowest_time){
                        lowest_time = time;
                    }
                }
            }
            curr_time = lowest_time;
            for(int i =0; i<end_times.size();i++){
                if(end_times[i] == lowest_time){
                    for(int j=0;j<4;j++){
                        curr_available_res[j] +=activities[list_repr[i]].resource_needed[j];
                    }
                    finish_act(activities[list_repr[i]]);
                }
            }
            step();
            return;
        }
        //Nemam dost zdroju, musim dokoncit nejblizsi a posunout cas
        for(int ii=0;ii<4;ii++){
            if(activities[act_to_plan].resource_needed[ii] > curr_available_res[ii]){

                int lowest_time = INT32_MAX;
                for(auto time : end_times){
                    if(time > curr_time){
                        if(time<lowest_time){
                            lowest_time = time;
                        }
                    }
                }
                curr_time = lowest_time;
                for(int i =0; i<end_times.size();i++){
                    if(end_times[i] == lowest_time){
                        for(int j=0;j<4;j++){
                            curr_available_res[j] +=activities[list_repr[i]].resource_needed[j];
                        }
                        finish_act(activities[list_repr[i]]);
                    }
                }
                step();
                return;
            }
        }
        // Constraints neporusim, planuji aktivitu
        schedule_activity(activities[act_to_plan]);
        
    }

    //Dokonceni aktivity, kontrola, ktere jsem uvolnil pro naplanovani.
    void finish_act(activity planned){
        
        currently_active.erase(planned.id);
        bool insert = true;
        for(int i=0; i<planned.num_succ;i++){
            for(auto pre : pred(activities[planned.succs[i]])){
                bool is_finished = true;
                for(auto curr_act : currently_active){
                    if(curr_act == pre){
                        is_finished = false;
                    }
                }
                
                bool was_planned = false;
                for(auto done : list_repr){
                    if(done == pre){
                        was_planned = true;
                    }
                }
                
                insert = was_planned && insert && is_finished;
            }
            if(insert){
                available_acts.insert(planned.succs[i]);
            }
            insert = true;
        }
    }

    //Naplanovani aktivity. Pokud je to prvni dummy, rovnou dokonci.
    void schedule_activity(activity planned){
        list_repr.emplace_back(planned.id);
        currently_active.insert(planned.id);
        for(int i=0;i<4;i++){
            curr_available_res[i] = curr_available_res[i] - planned.resource_needed[i];
        }
        
        times.emplace_back(curr_time);
        end_times.emplace_back(curr_time+planned.duration);
        available_acts.erase(planned.id);
        if(planned.duration == 0){
            finish_act(planned);
        }
        curr_index++;
        
    }

    //Vypocita seznam vsech predchudcu
    std::set<int> pred (activity a){
        std::set<int> pre;
        for(auto act : activities){
            if(act.id != a.id){
                for(int i =0; i<act.num_succ; i++){
                    if(act.succs[i] == a.id){
                        //std::cout << act.id << " " <<  std::endl;
                        pre.insert(act.id);
                        for(auto id : pred(act)){
                             pre.insert(id);
                        }
                    }
                }
            }
            
        }
        return pre;
    }
};

//Trida pro celou proceduru predstavujici algoritmus
class Procedure {
    public: 
        std::vector<int> spt_res;
        int spt_time;
        std::vector<int> best_schedule;
        int best_time;
        std::vector<int> curr_list;
        int curr_time;
        int N_0;
        int h;
        float T_0max;
        float alpha;
        int S;
        int C;
        float T;
        int N_S;
        int CP;
        int delta_max;
        int N_tot = 0;
        int N_avg;

        activity activities [62];  //Potreba zmenit pro zmenu datasetu J30 = [32], J60 = [62]
        //Hlavni funkce, implementujici predstaveny algoritmus
        void start(){
            init();
            
            for(int c=0; c<C;c++){
                T = T_0max;
                N_S = N_0;
                for(int s = 1; s<S+1;s++){
                    N_S = N_S * (1+h*s);
                    for(int n_p = 0; n_p<N_S; n_p++){
                        SA sa;
                        int i=0;
                        for(activity act : activities){
                            sa.activities[i] = act;
                            i++;
                        }
                        for(auto act : curr_list){
                            sa.given_list_repr.emplace_back(act);
                        }
                        for(int ind = 0; ind<4;ind++){
                            sa.resource_avaiability[ind] = resource_avaiabilities[ind];
                            sa.curr_available_res[ind] = resource_avaiabilities[ind];
                        }
                        
                        sa.compute_timespan();
                        N_tot++;
                        if(N_tot >= 20000){
                            return;
                        }
                        int delta =  sa.times[num_of_act-1]- curr_time;
                        if(delta<=0){
                            if(delta < 0){
                            }
                            curr_list.clear();
                            for(auto act : sa.given_list_repr){
                                curr_list.emplace_back(act);
                            }
                            curr_time = sa.times[num_of_act-1];
                            if(curr_time< best_time){
                                best_time = curr_time;
                                best_schedule.clear();
                                for(auto act : sa.given_list_repr){
                                    best_schedule.emplace_back(act);
                                }
                            }
                            if(best_time == CP){
                                return;
                            }
                        } else {
                            int random_num = rand() % 1000;
                            float y_rand = (float)random_num/1000;
                            if(pow(E, -(float)delta/T )> y_rand){
                                curr_list.clear();
                                for(auto act : sa.given_list_repr){
                                    curr_list.emplace_back(act);
                                }
                                curr_time = sa.times[num_of_act-1];
                            }
                        }
                    }
                    T = T*alpha;
                }
            }
            for(int n = 0; n<N_avg;n++){
                SA sa;
                int i=0;
                for(activity act : activities){
                    sa.activities[i] = act;
                    i++;
                }
                for(auto act : best_schedule){
                    sa.given_list_repr.emplace_back(act);
                }
                for(int ind = 0; ind<4;ind++){
                    sa.resource_avaiability[ind] = resource_avaiabilities[ind];
                    sa.curr_available_res[ind] = resource_avaiabilities[ind];
                }
                sa.compute_timespan();
                int delta =  sa.times[num_of_act-1]- best_time;
                if(delta <  0){
                    best_time =  sa.times[num_of_act-1];
                    best_schedule.clear();
                    for(auto act : sa.given_list_repr){
                        best_schedule.emplace_back(act);
                    }
                }
            }
            std::cout << "Ntot = " << N_tot << std::endl;   
        }
    private:
        //Inicializace parametru. Potreba upravit pro ruzne experimenty
        void init(){
            alpha = 0.25;
            S = 5;
            C = 1;
            h = 1;
            N_0=2;
            N_avg = 1000;
            best_time = spt_time;
            for(auto act : spt_res){
                best_schedule.emplace_back(act);
                curr_list.emplace_back(act);
            }
            curr_time = spt_time;
            delta_max = (float)spt_time/2;
            T_0max = 0.2*spt_time;
        }
};

//Trida pro vypocet inicialniho reseni dle SPT heuristiky
class SPT {
  public:             
    int resource_avaiability [4];
    activity activities [62];  //Potreba zmenit pro zmenu datasetu J30 = [32], J60 = [62]
    std::vector<int> list_repr;
    std::vector<int> times;
    std::vector<int> end_times;
    std::set<int> currently_active;
    std::set<int> available_acts;
    int curr_time;
    int curr_available_res[4];

    //Hlavni funkce
    void test (){
        init();
        while(list_repr.size()<num_of_act){
            step();
        }
    }

    //Vypocita CP value, coz je reseni, ktere ignoruje pozadavkz na zdroje.
    int compute_CP(){
        init();
        while(list_repr.size()<num_of_act){
            step_CP();
        }
        return list_repr[num_of_act-1];
    }


  private:  
    //Planovani a dokonceni prvni dummy activity a inicializace parametru.
    void init (){
        curr_time=0;
            
        for(int i = 0; i<4;i++){
            curr_available_res[i] = resource_avaiability[i];
        }
        available_acts.insert(0);
        schedule_activity(activities[0]);
    }

    //Krok simulace bez pozadavku na zdroje
    void step_CP(){
        int act_to_plan = pick_next();
        
        if(act_to_plan == -1){
            int lowest_time = INT32_MAX;
            for(auto time : end_times){
                if(time > curr_time){
                    if(time<lowest_time){
                        lowest_time = time;
                    }
                }
            }
            curr_time = lowest_time;
            for(int i =0; i<end_times.size();i++){
                if(end_times[i] == lowest_time){
                    finish_act(activities[list_repr[i]]);
                }
            }
            step_CP();
            return;
        } else {
            schedule_activity(activities[act_to_plan]);
        }
    }

    //Krok simulace
    void step(){
        //Vybere se dalsi aktivita k naplanovani
        int act_to_plan = pick_next();
        
        if(act_to_plan == -1){
            //Vsichni predchudci nejsou dokonceni, potreba posunout cas na nejblizsi konec aktivity
            int lowest_time = INT32_MAX;
            for(auto time : end_times){
                if(time > curr_time){
                    if(time<lowest_time){
                        lowest_time = time;
                    }
                }
            }
            curr_time = lowest_time;
            for(int i =0; i<end_times.size();i++){
                if(end_times[i] == lowest_time){
                    for(int j =0; j<4;j++){
                        curr_available_res[j] +=activities[list_repr[i]].resource_needed[j];
                    }
                    finish_act(activities[list_repr[i]]);
                }
            }
            step();
            return;
        }
        
        for(int ii=0; ii<4;ii++){
            if(activities[act_to_plan].resource_needed[ii] > curr_available_res[ii]){
                 //Nemam dost zdroju, potreba posunout cas na nejblizsi konec aktivity

                int lowest_time = INT32_MAX;
                for(auto time : end_times){
                    if(time > curr_time){
                        if(time<lowest_time){
                            lowest_time = time;
                        }
                    }
                }
                curr_time = lowest_time;
                for(int i =0; i<end_times.size();i++){
                    if(end_times[i] == lowest_time){
                        for(int j =0; j<4;j++){
                            curr_available_res[j] +=activities[list_repr[i]].resource_needed[j];
                        }
                        finish_act(activities[list_repr[i]]);
                    }
                }
                step();
                return;
            }
        }
        //Naplanovani aktivity
        schedule_activity(activities[act_to_plan]);
    }

    //Dokonceni aktivity, kontrola, ktere jsem dokoncenim uvolnil
    void finish_act(activity planned){
        currently_active.erase(planned.id);
        bool insert = true;
        for(int i=0; i<planned.num_succ;i++){
            for(auto pre : pred(activities[planned.succs[i]])){
                bool is_finished = true;
                for(auto curr_act : currently_active){
                    if(curr_act == pre){
                        is_finished = false;
                    }
                }
                
                bool was_planned = false;
                for(auto done : list_repr){
                    if(done == pre){
                        was_planned = true;
                    }
                }
                
                insert = was_planned && insert && is_finished;
            }
            if(insert){
                available_acts.insert(planned.succs[i]);
            }
            insert = true;
        }
    }
    //Vybrani dalsi aktivity k naplanovani dle SPT
    int pick_next(){
        
        int res_id;
        int lowest_time;
        int i =0;
        if(available_acts.empty()){
            return -1;
        }
        for(auto act : available_acts){
            if(i == 0){
                lowest_time = activities[act].duration;
                res_id = act;
            }
            if(activities[act].duration < lowest_time){
                lowest_time = activities[act].duration;
                res_id = act;
            }
            i++;
        }
        return res_id;
    }
    //Planovani aktivity
    void schedule_activity(activity planned){
        list_repr.emplace_back(planned.id);
        currently_active.insert(planned.id);
        for(int i =0; i<4;i++){
            curr_available_res[i] = curr_available_res[i] - planned.resource_needed[i];
        }
        times.emplace_back(curr_time);
        end_times.emplace_back(curr_time+planned.duration);
        available_acts.erase(planned.id);
        if(planned.duration == 0){
            finish_act(planned);
        }
        
    }

    //Vypocita sadu predchudcu.
    std::set<int> pred (activity a){
        std::set<int> pre;
        for(auto act : activities){
            if(act.id != a.id){
                for(int i =0; i<act.num_succ; i++){
                    if(act.succs[i] == a.id){
                        pre.insert(act.id);
                        for(auto id : pred(act)){
                             pre.insert(id);
                        }
                    }
                }
            }
        }
        return pre;
    }
};



//Debug funkce
void print_act(activity a){
    std::cout << "id: " << a.id << std::endl; 
    std::cout << "duration: " << a.duration << std::endl; 
    std::cout << "num succs: " << a.num_succ << std::endl; 
    int i=0;
    for(auto rec_ned : a.resource_needed) {
        std::cout << "rec_needed" << i << ": " << rec_ned << std::endl; 
        i++;
    }
    for (i =0; i< a.num_succ; i++){
        std::cout << "succ" << i << ": " << a.succs[i] << std::endl; 
    }
}

int main() {
    int opt_res = 0;
    float var=0;
    int worse = 0;
    int total_num = 0;
    srand (time(NULL));
    //Nacitani optimalniho reseni
    std::ifstream hrsfile("j60hrs.sm"); //Potreba zmenit pro zmenu datasetu J30 = j30hrs.sm, J60 = j60hrs.sm
    int opt[480];
    if (hrsfile.is_open()) {
        std::string line;
        int i =0;
        while (std::getline(hrsfile, line)) {
            if(i>3){
                int numofnums=0;
                int index1;
                int index2;
                int makespan;
                std::string num = "";
                for(auto a : line){
                    
                    if(a=='	' && numofnums==0){
                        index1 = stoi(num);
                        numofnums++;
                        num = "";
                    } else if (a=='	' && numofnums==1){
                        index2 = stoi(num);
                        numofnums++;
                        num = "";
                    } else if(a=='	' && numofnums==2){
                        makespan = stoi(num);
                        numofnums++;
                        num = "";
                    } else {
                        num+=a;
                    }
                }
                opt[(index1-1)*10 + index2-1] = makespan;
            }
            i++;
        }
    }
    std::string path = "j60rcp"; //Potreba zmenit pro zmenu datasetu J30 = j30rsp, J60 = j60rsp
    int num = 62;  //Potreba zmenit pro zmenu datasetu J30 = 32, J60=62
    //Nacitani datasetu a vypocet
    for (const auto & entry : directory_iterator(path)){
        int i = 0;
        activity activities[num];

        std::cout << entry.path() << std::endl;
        int ch_index = 0;
        int index1;
        int index2;
        std::string pa =  entry.path();
        std::string num = "";
        for(auto ch: pa){
            if(ch_index > 9){
                if(ch == '_'){                  
                    index1 = stoi(num);
                    num = "";
                } else if(ch == '.'){
                    index2 = stoi(num);
                    break;
                } else {
                    num+=ch;
                }
            }
            ch_index++;
        }
        int fin_index = (index1-1)*10 + index2-1;
        std::ifstream file(entry.path());
        if (file.is_open()) {
            std::string line;
            while (std::getline(file, line)) {
                if (i ==0) {
                    std::string num = "";
                    for (auto a : line){
                        if(a == ' '){
                            num_of_act = stoi(num);
                            break;
                        } else {
                            num += a;
                        }
                    }
                } else if (i == 1){
                    int ii = 0;
                    std::string num = "";
                    for (auto a : line){
                        if(a == ' '){
                            if(num != ""){
                                resource_avaiabilities[ii] = stoi(num);
                                num = "";
                                ii++;
                            }

                        } else {
                            num += a;
                        }
                    }
                } else {
                    activities[i-2].id = i-2;
                    std::string num = "";
                    int ii = 0;
                    for (auto a : line){
                        if(a == ' '){
                            if(num != ""){
                                if(ii == 0){
                                    activities[i-2].duration = stoi(num);
                                } else if(ii ==1){
                                    activities[i-2].resource_needed[0] = stoi(num);
                                } else if(ii ==2){
                                    activities[i-2].resource_needed[1] = stoi(num);
                                } else if(ii ==3){
                                    activities[i-2].resource_needed[2] = stoi(num);
                                } else if(ii ==4){
                                    activities[i-2].resource_needed[3] = stoi(num);
                                } else if(ii ==5){
                                    activities[i-2].num_succ = stoi(num);
                                } else {
                                    activities[i-2].succs[ii-6] = stoi(num)-1;
                                }
                                num = "";
                                ii++;
                            }

                        } else {
                            num += a;
                        }
                    }
                }
                i++;
            }
            file.close();
        }
        
        SPT spt;
        i=0;
        for(activity act : activities){
            spt.activities[i] = act;
            i++;
        }
        
        spt.curr_time = 0;
        for(int j=0;j<4;j++){
            spt.curr_available_res[j] = resource_avaiabilities[j];
            spt.resource_avaiability[j] = resource_avaiabilities[j]; 
        }
        spt.test();
        std::cout<< "SPT time = " << spt.times[num_of_act-1] << std::endl;
        SPT cp_spt;
        i=0;
        for(activity act : activities){
            cp_spt.activities[i] = act;
            i++;
        }
        
        cp_spt.curr_time = 0;
        for(int j=0;j<4;j++){
            cp_spt.curr_available_res[j] = resource_avaiabilities[j];
            cp_spt.resource_avaiability[j] = resource_avaiabilities[j];
        }
        int cp_val = cp_spt.compute_CP();
        SA sa;

        i=0;
        for(activity act : activities){
            sa.activities[i] = act;
            i++;
        }

        i=0;
        for(auto act : spt.list_repr){
            sa.given_list_repr.emplace_back(act);
            i++;
        }
        for(int j=0;j<4;j++){
            sa.resource_avaiability[j] = resource_avaiabilities[j];
            sa.curr_available_res[j] = resource_avaiabilities[j];
        }
        
        sa.compute_timespan();
        

        Procedure procedure;
        for(auto act : spt.list_repr){
            procedure.spt_res.emplace_back(act);
            i++;
        }
        i=0;
        for(activity act : activities){
            procedure.activities[i] = act;
            i++;
        }
        procedure.CP = cp_val;
        procedure.spt_time = spt.times[num_of_act-1];
        procedure.start();
        std::cout << "best time = " << procedure.best_time << " optimal = " << opt[fin_index]<<std::endl;
        if(opt[fin_index] != procedure.best_time){
            if(procedure.best_time<opt[fin_index]){
                worse++;
            }
            var += abs(procedure.best_time - opt[fin_index]);
        } else {
            opt_res++;
        }
        total_num++;
        std::cout<< "total num " << total_num << std::endl;
        std::cout<< ",variance = " << var << "     ";
        std::cout<< ",opt count = " << opt_res << std::endl;
    }
    std::cout<< "better = " << worse << std::endl;
    std::cout<< "fin variance = " << var/480 << std::endl;
	return 0;
}