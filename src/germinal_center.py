import random
from src.constants import *
from scipy.spatial.distance import hamming
from deap import tools

class GerminalCenter():
    def __init__(self, toolbox, antigen, antibodies):
        self.toolbox = toolbox
        self.antigen = antigen
        self.antibodies = antibodies
        self._init_antibodies()
        self.epitopes = self._gen_epitopes(EPITOPES_PER_ANTIGEN)

    def _init_antibodies(self):
        for i in range(len(self.antibodies)):
            self.antibodies[i].fitness.values = self.toolbox.init_antibody_fitness(self.antibodies[i])

    def _bind(self, antibody, epitope):
        #FIXME do better when they have different lengths
        if len(antibody) != len(epitope):
            return False, 9999999
        dist = hamming(antibody, epitope)
        if dist < HAM_DIST:
            return True, dist
        return False, dist


    def evolve(self):
        #Each antibody tries to bind to an epitope
        some_did_not_bind = False
        num_bound = 0
        for antibody in self.antibodies:
            #Try to bind
            epitope = random.choice(self.epitopes)
            did_bind, dist = self._bind(antibody, epitope)
            if did_bind:
                fit = antibody.fitness.values
                #new_fit = (fit[0]+1,)
                new_fit = (dist,)
                antibody.fitness.values = new_fit
                print(antibody.fitness.values)
                num_bound += 1
            else:
                some_did_not_bind = True

        if not some_did_not_bind:
            print("All binded!")
            return True
        if num_bound > 0:
            print(f"{num_bound} bounded")

        next_gen = []

        #Pick the top X performing antibodies to clone/mutate
        num_to_keep = int(ANTIBODY_CLONE_PCT * NUM_ANTIBODIES)
        bests = tools.selBest(self.antibodies, k=num_to_keep)
        next_gen = []
        for best in bests:
            num_clones_per = int(NUM_ANTIBODIES / num_to_keep) - 1
            next_gen.append(best)
            for clone_num in range(num_clones_per):
                clone = self.toolbox.clone(best)
                self._mutate(clone)
                next_gen.append(clone)
        assert len(next_gen) == NUM_ANTIBODIES, f"expected {NUM_ANTIBODIES} antibodies, had {len(next_gen)}"
        random.shuffle(next_gen)
        self.antibodies = next_gen
        return False

    def _mutate(self, clone):
        #Flip a random bit, modify in-place
        ind = random.randint(0, len(clone) - 1)
        if clone[ind] == 0:
            clone[ind] = 1
        else:
            clone[ind] = 0
    
    def _gen_epitopes(self, num_epitopes):
        epitopes = []
        for _ in range(num_epitopes):
            length_choices = [36, 42, 48]
            epitope_len = random.choice(length_choices)
            start_ind = random.randint(0, len(self.antigen) - epitope_len)
            epitope = self.antigen[start_ind : start_ind + epitope_len]
            epitopes.append(epitope)
        return epitopes