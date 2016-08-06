class PolynomialKontsevichAutomorphism(SageObject):
    def __init__(self):

        data_dir = "/Users/dylanrupel/.sage/bin/data/"
        computation_dir = "/Users/dylanrupel/Documents/Math/Research/Computations/"

        S.<v> = ZZ[]
        self.S = S
        self.Q = S.fraction_field()
        base_ring.<j11,j21,j22,j23,j24> = Q[]
        self.base_ring = base_ring

        """
        Change base ring to correspond to changes here
        """
        self.P1 = [[0,1],[1,1]] #formal polynomial, first coordinate is the exponent, second coordinate is the coefficient
        self.P1rev = [[0,1],[1,1]]
        self.P2 = [[0,1],[1,j23],[2,j22],[3,j21],[4,1]]
        self.P2rev = [[0,1],[1,j21],[2,j22],[3,j23],[4,1]]

        self.comm_matrix = matrix([[0,1],[-1,0]])
        self.comm_matrix.set_immutable()
        self.A = BasedQuantumTorus(self.base_ring,2,self.Q(v),self.comm_matrix)
        self.B = self.A.basis()
        self.X1 = self.B[(1,0)]
        self.X1i = self.B[(-1,0)]
        self.X2 = self.B[(0,1)]
        self.X2i = self.B[(0,-1)]
        
        F.<X,Y,Z> = FreeAlgebra(self.base_ring,4)
        self.X = X
        self.Y = Y
        self.Z = Z
        self.F = F
        self.M = self.F.monoid()
        self.Xi = self.M(1)
        self.Yi = self.M(1)
        self.Zi = self.M(1)
        self.Xi._element_list.append((1,-1))
        self.Yi._element_list.append((2,-1))
        self.Zi._element_list.append((3,-1))
        self.Xi = self.F(self.Xi)
        self.Yi = self.F(self.Yi)
        self.Zi = self.F(self.Zi)
        
        self._X = matrix(QQ,[[1,3],[0,1]])
        self._Xi = self._X.inverse()
        self._Y = matrix(QQ,[[1,2],[0,1]])
        self._Yi = self._Y.inverse()
        self._Z = matrix(QQ,[[0,0],[0,0]])
        for i,j in self.P1:
            self._Z += j*self._Y^i
        self._Zi = self._Z.inverse()
    
    def quantum_multinomial(self, n,bfk,q):
        """
        Input:
            -Integer n
            -List of integers bfk
            -An element q of some ring
        """
        sum = 0
        for k in bfk:
            sum += k
            if k > n:
                return 0
            if k < 0:
                return 0
        if sum > n:
            return 0
        if sum == 0:
            return 1
        output = 0
        for i in range(0,bfk.__len__()+1):
            exp = 0
            for j in range(i+1,bfk.__len__()+1):
                exp += (j-i)*bfk[j-1]
            if i == 0:
                output += q^exp*quantum_multinomial(n-1,bfk,q)
            else:
                bfk[i-1] -= 1
                output += q^exp*quantum_multinomial(n-1,bfk,q)
                bfk[i-1] += 1
        return output
    

    def latex_noncomm_elt(self, elt, Factored=true):
        output = ""

    
    def compute_gen_noncomm_variables(self):
    
            working_dir = "Generalized_Cluster_Variables/NonCommutative/"
            filename = "d1="+str(self.P1[self.P1.__len__()-1][0])+",d2="+str(self.P2[self.P2.__len__()-1][0])+".tex"
            TeXFile=open(computation_dir+working_dir+filename,'w')
            TeXFile.write("\\documentclass{article}\n")
            TeXFile.write("\\usepackage{amsmath, amssymb, latexsym}\n")
            TeXFile.write("\\usepackage[margin=1in]{geometry}\n\n")
            TeXFile.write("\\begin{document}\n\n")
            TeXFile.write("Let $P_1=")
            self.P1_str = ""
            for term in self.P1:
                if term[0] == 0:
                    self.P1_str += "1"
                else:
                    if term[1] > 1:
                        self.P1_str += "+"+str(term[1])
                    elif term[1] == 1:
                        self.P1_str += "+"
                    elif self.P1_str < 0:
                        self.P1_str += str(term[1])
                    if term[1] != 0:
                        self.P1_str += "z"
                        if term[0] > 1:
                            self.P1_str += "^{"+str(term[0])+"}"
            TeXFile.write(self.P1_str+"$ and $P_2=")
            self.P2_str = ""
            for term in self.P2:
                if term[0] == 0:
                    self.P2_str += "1"
                else:
                    if term[1] > 1:
                        self.P2_str += "+"+str(term[1])
                    elif term[1] == 1:
                        self.P2_str += "+"
                    elif self.P2_str < 0:
                        self.P2_str += str(term[1])
                    if term[1] != 0:
                        self.P2_str += "z"
                        if term[0] > 1:
                            self.P2_str += "^{"+str(term[0])+"}"
            TeXFile.write(self.P2_str+"$.\\\\\n\n")
            for n in range(0,7):
                TeXFile.write("$Y_{"+str(n+1)+"}=")
                noncomm_var = self.fsPKont(self.P1,self.P2,n)
                TeXFile.write(str(noncomm_var)+"$\\\\\n\n")
    
            TeXFile.write("\\end{document}")
            TeXFile.close()
            
            import subprocess
            subprocess.call(['pdflatex', '-halt-on-error', filename], cwd=computation_dir+working_dir, stdout=subprocess.PIPE)
    
    
    def compute_gen_quantum_variables(self):
    
            working_dir = "Generalized_Cluster_Variables/Quantum/"
            filename = "d1="+str(self.P1[self.P1.__len__()-1][0])+",d2="+str(self.P2[self.P2.__len__()-1][0])+".tex"
            TeXFile=open(computation_dir+working_dir+filename,'w')
            TeXFile.write("\\documentclass{article}\n")
            TeXFile.write("\\usepackage{amsmath, amssymb, latexsym}\n")
            TeXFile.write("\\usepackage[margin=1in]{geometry}\n\n")
            TeXFile.write("\\begin{document}\n\n")
            TeXFile.write("Let $P_1=")
            self.P1_str = ""
            for term in self.P1:
                if term[0] == 0:
                    self.P1_str += "1"
                else:
                    if term[1] > 1:
                        self.P1_str += "+"+str(term[1])
                    elif term[1] == 1:
                        self.P1_str += "+"
                    elif self.P1_str < 0:
                        self.P1_str += str(term[1])
                    if term[1] != 0:
                        self.P1_str += "z"
                        if term[0] > 1:
                            self.P1_str += "^{"+str(term[0])+"}"
            TeXFile.write(self.P1_str+"$ and $P_2=")
            self.P2_str = ""
            for term in self.P2:
                if term[0] == 0:
                    self.P2_str += "1"
                else:
                    if term[1] > 1:
                        self.P2_str += "+"+str(term[1])
                    elif term[1] == 1:
                        self.P2_str += "+"
                    elif self.P2_str < 0:
                        self.P2_str += str(term[1])
                    if term[1] != 0:
                        self.P2_str += "z"
                        if term[0] > 1:
                            self.P2_str += "^{"+str(term[0])+"}"
            TeXFile.write(self.P2_str+"$.\\\\\n\n")
            for n in range(0,7):
                TeXFile.write("$X_{"+str(n+1)+"}=")
                noncomm_var = self.fsPKont(self.P1,self.P2,n)
                quantum_var = quantum_specialization(noncomm_var,self.P1)
                TeXFile.write(self.A.latex_element(quantum_var)+"$\\\\\n\n")
    
            TeXFile.write("\\end{document}")
            TeXFile.close()
            
            import subprocess
            subprocess.call(['pdflatex', '-halt-on-error', filename], cwd=computation_dir+working_dir, stdout=subprocess.PIPE)
    
    
    def compute_gen_variables_coeff(self):
    
            working_dir = "Generalized_Cluster_Variables/Quantum/"
            d1 = self.P1[self.P1.__len__()-1][0]
            d2 = self.P2[self.P2.__len__()-1][0]
            filename = "d1="+str(d1)+",d2="+str(d2)+"-coefficients.tex"
            TeXFile=open(computation_dir+working_dir+filename,'w')
            TeXFile.write("\\documentclass{article}\n")
            TeXFile.write("\\usepackage{amsmath, amssymb, latexsym}\n")
            TeXFile.write("\\usepackage[margin=1in]{geometry}\n\n")
            TeXFile.write("\\newcommand{\\JJ}{\\mathbb{J}}\n\n")
            TeXFile.write("\\begin{document}\n\n")
            TeXFile.write("Let $P_1=")
            self.P1_str = ""
            for term in self.P1:
                if term[0] == 0:
                    self.P1_str += "1"
                else:
                    if term[1] != 0:
                        self.P1_str += "+z"
                        if term[0] > 1:
                            self.P1_str += "^{"+str(term[0])+"}"
            TeXFile.write(self.P1_str+"$ and $P_2=")
            self.P2_str = ""
            for term in self.P2:
                if term[0] == 0:
                    self.P2_str += "1"
                else:
                    if term[1] != 0:
                        self.P2_str += "+z"
                        if term[0] > 1:
                            self.P2_str += "^{"+str(term[0])+"}"
            TeXFile.write(self.P2_str+"$.\\\\\n\n")
            for n in range(0,7):
                TeXFile.write("$X_{"+str(n+1)+"}=")
                noncomm_var = self.fsPKont(self.P1,self.P2,n)
                quantum_var = quantum_specialization(noncomm_var,self.P1)
    
                sup = quantum_var.support()
                rank_vector = [-sup[0][0],-sup[0][1]]
                sup.reverse()
                rep_coeffs = quantum_var.monomial_coefficients()
                output = ""
                for key in sup:
                    if output!="":
                        output = output+"+"
                
                    rep_coeff = rep_coeffs[key]
                    rep_dict = rep_coeff._mpoly_dict_recursive()
                    rep_dict_keys = rep_dict.keys()
                    if rep_dict_keys.__len__() > 1:
                        output += "\\big["
                    sum_switch = False
                    for rep_key in rep_dict_keys:
                        if sum_switch:
                            output += "+"
                        vcoeff = rep_dict[rep_key]
                        if vcoeff != self.base_ring.one():
                            output += self.A.latex_coeff(vcoeff)
                        for i in range(0,d1):
                            if rep_key[i] > 0:
                                output += "\\JJ_{1,"+str(i+1)+"}"
                                if rep_key[i] > 1:
                                    output += "^{"+str(rep_key[i])+"}"
                        for i in range(0,d2):
                            if rep_key[d1+i] > 0:
                                output += "\\JJ_{2,"+str(i+1)+"}"
                                if rep_key[d1+i] > 1:
                                    output += "^{"+str(rep_key[d1+i])+"}"
                        sum_switch = True
                    if rep_dict_keys.__len__() > 1:
                        output += "\\big]"
                    output += "X^{("
                    for i in range(0,2):
                        if i != 0:
                            output += ","
                        output += str(key[i])
                    output += ")}"
                TeXFile.write(output+"$\\\\\n\n")
    
            TeXFile.write("\\end{document}")
            TeXFile.close()
            
            import subprocess
            subprocess.call(['pdflatex', '-halt-on-error', filename], cwd=computation_dir+working_dir, stdout=subprocess.PIPE)
    
    
    def quantum_specialization(self, element,pol):
        """
        Input:
            -an element of F
            -a formal polynomial self.P1
        """
            
        if element.parent() == self.M:
            output = self.A(1)
            for (i,j) in element._element_list:
                if i == 1:
                    if j > 0:
                        output = output*(v*self.X1)^j
                    if j < 0:
                        output = output*(v^(-1)*self.X1i)^(-j)
                if i == 2:
                    if j > 0:
                        output = output*(v^(-1)*self.X2)^j
                    if j < 0:
                        output = output*(v*self.X2i)^(-j)
                if i == 3:
                    if j > 0:
                        factor = self.A(0)
                        for term in pol:
                            if term[0] != 0:
                                factor += term[1]*(v^(-1)*self.X2)^term[0]
                            else:
                                factor += self.base_ring.gens()[pol.__len__() - 2]
                        output = output*factor^j
            return output
            
        if element.parent() == self.F:
            output = self.A(0)
            if element == self.A(0):
                return output
            mon_coeff = element._monomial_coefficients
            for key in mon_coeff.keys():
                output = output + mon_coeff[key]*quantum_specialization(key,pol)
    
    
            d1 = self.P1[self.P1.__len__()-1][0]
            d2 = self.P2[self.P2.__len__()-1][0]
            sup = output.support()
            rank_vector = [-sup[0][0],-sup[0][1]]
            if rank_vector[0] < 0 or rank_vector[1] < 0:
                return v^(-1)*output
            else:
                sup.reverse()
                rep_coeffs = output.monomial_coefficients()
                proper_output = self.A(0)
                for key in sup:
                    rep_coeff = rep_coeffs[key]
                    rep_dict = rep_coeff._mpoly_dict_recursive()
                    rep_dict_keys = rep_dict.keys()
                    for rep_key in rep_dict_keys:
                        term = rep_dict[rep_key]
                        dim1_sum = 0
                        for i in range(0,d1):
                            dim1_sum += rep_key[i]*(i+1)
                            term = term*self.base_ring.gens()[i]^rep_key[i]
                        dim = d1*rank_vector[0] - (key[1] + rank_vector[1]) - dim1_sum
                        if dim%d1 != 0:
                            print "rank=", rank_vector
                            print "term=", key
                            raise ValueError, "The difference "+str(dim)+" in dimensions is not divisible by "+str(d1)+"."
                        exp1_offset = dim/d1
                        if dim < 0:
                            raise ValueError, "The exponent "+str(exp1_offset)+" should be positive."
                        dim2_sum = 0
                        for i in range(0,d2):
                            dim2_sum += rep_key[d1+i]*(i+1)
                            term = term*self.base_ring.gens()[d1+i]^rep_key[d1+i]
                        dim = key[0] + rank_vector[0] - dim2_sum
                        if dim%d2 != 0:
                            print "rank=", rank_vector
                            print "key=", key
                            print "reps=", rep_key
                            raise ValueError, "The difference "+str(dim)+" in dimensions is not divisible by "+str(d2)+"."
                        exp2_offset = dim/d2
                        if dim < 0:
                            raise ValueError, "The exponent "+str(exp2_offset)+" should be positive."
                        proper_output += self.base_ring.gens()[d1 - 1]^exp1_offset*self.base_ring.gens()[d1+d2-1]^exp2_offset*term*self.B[key]
                return v^(-1)*proper_output
    
    
    def Pexternal_Kont(self, Kont_pair, P):
        temp = copy(Kont_pair)
        Kont_pair[0] = temp[0]*temp[1]*temp[0].inverse()
        factor = 0
        for i,j in P:
            factor += j*temp[1]^i
        Kont_pair[1] = factor*temp[0].inverse()
        return Kont_pair
    
    #use this for matrix computations
    def nPKont(self, P1, P2, n, _X, _Y):
        """
        #numerical Kontsevich automorphism
        Input: 
            -formal polynomials self.P1 and self.P2
            -positive integer n
            -invertible matrices _x, _y
        output:
            -a new pair of matrices Kont_pair
        """
        Kont_pair = [_X,_Y]
        for i in range(0,n):
            if i%2 == 0:
                Kont_pair = Pexternal_Kont(Kont_pair,P1)
            else:
                Kont_pair = Pexternal_Kont(Kont_pair,P2)
        return Kont_pair
    
    
    def fsPKont(self, n):
        #use this for factored symbolic computations,  z=(1+y^r) will not be expanded, relatively fast
        output = self.Y
        for i in range(0,n):
            if (n-i)%4 == 0:
                output = self.fsPKont_term(output,self.P2rev,self.P1)
            elif (n-i)%4 == 1:
                output = self.fsPKont_term(output,self.P1,self.P2)
            elif (n-i)%4 == 2:
                output = self.fsPKont_term(output,self.P2,self.P1rev)
            else:
                output = self.fsPKont_term(output,self.P1rev,self.P2rev)
        return output
    
    #use this for expanded symbolic computations, much slower
    def sPKont(self, n):
        output = self.Y
        for i in range(0,n):
            if (n-i)%4 == 0:
                output = self.sPKont_term(output,self.P2rev)
            elif (n-i)%4 == 1:
                output = self.sPKont_term(output,self.P1)
            if (n-i)%4 == 2:
                output = self.sPKont_term(output,self.P2)
            else:
                output = self.sPKont_term(output,self.P1rev)
        return output
    
    def fsPKont_term(self, element, P1, P2):
        """
        Input:
            -An element of F
            -A formal polynomial P1 used in this iteration
            -A formal polynomial P2 describing the Y-expansion of Z from the previous iteration
        """
        kont_X = self.X*self.Y*self.Xi
        kont_Xi = self.X*self.Yi*self.Xi
        kont_Y = self.Z*self.Xi
            
        if element.parent() == self.M:
            output = self.F(1)
            for (i,j) in element._element_list:
                if i == 1:
                    if j > 0:
                        output = output*(kont_X^j)
                    if j < 0:
                        output = output*(kont_Xi^(-j))
                if i == 2:
                    if j > 0:
                        output = output*(kont_Y^j)
                    if j < 0:
                        for k in range(0,-j):
                            output = output*self.X*self.Zi
                if i == 3:
                    if j > 0:
                        self.P2_factor = self.F(0)
                        for term in P2:
                            P2_factor += term[1]*kont_Y^term[0]
                        output = output*P2_factor^j
            return output
            
        if element.parent() == self.F:
            output = self.F(0)
            if element == self.F(0):
                return output
            mon_coeff = element._monomial_coefficients
            for key in mon_coeff.keys():
                output = output + mon_coeff[key]*self.fsPKont_term(key,P1,P2)
            
            return self.Psimplify(output,P1)
    
    def Pexpand(self, element, P):
        """
        Input:
            -An element of F
            -A formal polynomial P
        Output:
            -the element with all z's expanded
        """
    
        if element.parent() == self.M:
            output = self.F(1)
            for (i,j) in element._element_list:
                if i == 1:
                    if j > 0:
                        output = output*self.X^j
                    if j < 0:
                        output = output*self.Xi^(-j)
                if i == 2:
                    if j > 0:
                        output = output*self.Y^j
                    if j < 0:
                        output = output*self.Yi^(-j)
                if i == 3:
                    if j > 0:
                        P_factor = self.F(0)
                        for term in P:
                            P_factor += term[1]*self.Y^term[0]
                        output = output*P_factor^j
            return output
            
        if element.parent() == self.F:
            output = self.F(0)
            if element == self.F(0):
                return output
            mon_coeff = element._monomial_coefficients
            for key in mon_coeff.keys():
                output = output + mon_coeff[key]*Pexpand(key,P)
    
        return output
    
    
    def sPKont_term(self, element, P):
        """
        Input:
            -An element of F
            -A formal polynomial P
        Output:
            -the symbolic polynomial Kontsevich automorphism F_P applied to element
        """
        kont_X = self.X*self.Y*self.Xi
        kont_Xi = self.X*self.Yi*self.Xi
        kont_Y = self.F(0)
        for term in P:
            kont_Y += term[1]*self.Y^term[0]
        kont_Y = kont_Y*self.Xi
            
        if element.parent() == self.M:
            output = self.F(1)
            if self.Y_degree(element) < 0:
                raise ValueError, "The monomial %s does not have positive Y-degree."%element
            for (i,j) in element._element_list:
                if i == 1:
                    if j > 0:
                        output = output*(kont_X^j)
                    if j < 0:
                        output = output*(kont_Xi^(-j))
                if i == 2:
                    if j > 0:
                        output = output*(kont_Y^j)
                    if j < 0:
                        for k in range(0,-j):
                            output = output*self.X
                            #print "Kont_element",Kont_element
                            output = self.Polynomial_right_divide(output, P)
            return output
            
        if element.parent() == self.F:
            mon_coeff = element._monomial_coefficients
            
            #Find maximal non-negative Y prefixes    
            suffix_dict = dict()
            for key in mon_coeff.keys():
                div = self.Y_split(key)
                if suffix_dict.has_key(div[1]):
                    suffix_dict[div[1]] = suffix_dict[div[1]] + mon_coeff[key]*self.F(div[0])
                else:
                    suffix_dict.setdefault(div[1],mon_coeff[key]*self.F(div[0]))
                    
            #print "Applying Kont to non-negative Y degree prefixes"
                   
            Kont_dict = dict()
            for Key in suffix_dict.keys():
                entry = self.F(0)
                mon_coeff = suffix_dict[Key]._monomial_coefficients
                for key in mon_coeff.keys():
                    entry = entry + mon_coeff[key]*self.sPKont_term(key,P)
                Kont_dict.setdefault(Key,entry)
            
            #print "Kont applied to non-negative Y degree prefixes, moving to suffixes"
                
            Keys = Kont_dict.keys()
            Keys.sort()
            while Keys != [self.M(1)]:
                key = Keys.pop()
                entry = Kont_dict.pop(key)
                elt_list = key._element_list
                (i,j) = elt_list.pop(0)
                if i == 1:
                    if j > 0:
                        entry = entry*(kont_X^j)
                    if j < 0:
                        entry = entry*(kont_Xi^(-j))
                if i == 2:
                    if j > 0:
                        entry = entry*(kont_Y^j)
                    if j < 0:
                        entry = entry*self.X
                        entry = self.Polynomial_right_divide(entry, P)
                        elt_list.reverse()
                        elt_list.append((i,j+1))
                        elt_list.reverse()
                if Kont_dict.has_key(key):
                    Kont_dict[key] = Kont_dict[key] + entry
                else:
                    Kont_dict.setdefault(key,entry)
                    Keys.append(key)
                    Keys.sort()
            
            return Kont_dict[self.M(1)]
    
    def Y_degree(self, mon):
        deg = 0
        for (i,j) in mon._element_list:
            if i == 2:
                deg += j
        return deg
        
    def Z_split(self, mon):
        # splits the monomial mon into its maximal non-negative Z prefix and the corresponding suffix which begins with Z^{-1}
        prefix = self.M(1)
        suffix = self.M(1)
        switch = true
        for elt in mon._element_list:
            if elt[0] != 3 and switch:
                prefix._element_list.append(elt)
            else:
                if elt[1] >= 0 and switch:
                    prefix._element_list.append(elt)
                else:
                    switch = false
                    suffix._element_list.append(elt)
                    
        return [prefix,suffix]
        
    def commute(self, element):
        #element belongs to the monoid M, commute Y's and Z's, move Y's to the left
        if element == self.M(1):
            #element is the identity
            return element
        working_list = copy(element._element_list)
        output = self.M(1)
        while working_list != []:
            factor = working_list.pop(0)
            if factor[0] < 3 or working_list == []:
                #the factor is a power of X or Y
                output._element_list.append(factor)
            else:
                #the factor is a power of Z
                if working_list[0][0] == 1:
                    #the next factor is a power of X
                    output._element_list.append(factor)
                else:
                    #the next factor is a power of Y, so we should commute it past the Z
                    output._element_list.append(working_list.pop(0))
                    if working_list != []:
                        #combine with the next term if possible
                        if working_list[0][0] == 1:
                            #the next factor is a power of X
                            output._element_list.append(factor)
                        else:
                            #the next factor is a power of Z
                            new_term = (factor[0],factor[1]+working_list[0][1])
                            working_list.pop(0)
                            if new_term[1] != 0:
                                output._element_list.append(new_term)
                    else:
                        #no combining necessary
                        output._element_list.append(factor)
        return output
        
    def Psimplify(self, element,P):
        """
        Input:
            -an element of F
            -a formal polynomial P describing the Y-expansion of Z
        """
        if element.parent() == self.M:
            return self.commute(element)
                
        if element.parent() == self.F:
            #element is contained in the free algebra F, simplify each monomial then remove Z^{-1}
            output = self.F(0)
            if element == self.F(0):
                return output
            mon_coeff = element._monomial_coefficients
            
            #simplify each monomial then find maximal non-negative z prefixes    
            suffix_dict = dict()
            for key in mon_coeff.keys():
                simp_key = self.commute(key)
                div = self.Z_split(simp_key)
                if suffix_dict.has_key(div[1]):
                    suffix_dict[div[1]] = suffix_dict[div[1]] + mon_coeff[key]*self.F(div[0])
                else:
                    suffix_dict.setdefault(div[1],mon_coeff[key]*self.F(div[0]))
                
            Keys = suffix_dict.keys()
            Keys.sort()
            while Keys != [self.(1)]:
                key = Keys.pop()
                entry = suffix_dict.pop(key)
                elt_list = key._element_list
                if elt_list[0][0] == 3 and elt_list[0][1] < 0:
                    (i,j) = elt_list.pop(0)
                    entry = self.Polynomial_right_divide(entry, P)
                    if j+1 != 0:
                        elt_list.reverse()
                        elt_list.append((i,j+1))
                        elt_list.reverse()
                else:
                    div = self.Z_split(key)
                    entry = entry*self.F(div[0])
                    key = div[1]
                if suffix_dict.has_key(key):
                    suffix_dict[key] += entry
                else:
                    suffix_dict.setdefault(key,entry)
                    Keys.append(key)
                    Keys.sort()
            
            return suffix_dict[self.M(1)]
    
    def Polynomial_right_divide(self, element,P):
        """
        Input:
            -an element of F
            -a formal polynomial P describing the Y-expansion of Z
        """
        output = self.F(0)
        if element == self.F(0):
            return output
        mon_coeff = element._monomial_coefficients
        #print "mon_coeff",mon_coeff
        trailing_Y_dict = dict()
        for key in mon_coeff.keys():
            elt_list = key._element_list
            #print "elt_list",elt_list
            if elt_list == [] or elt_list[len(elt_list)-1][0] != 2:
                Y_exp = 0
            else:
                Y_exp = elt_list[len(elt_list)-1][1]
            if trailing_Y_dict.has_key(Y_exp):
                trailing_Y_dict[y_exp] = trailing_Y_dict[Y_exp] + mon_coeff[key]*self.F(key)
            else:
                trailing_Y_dict.setdefault(Y_exp,mon_coeff[key]*self.F(key))
    
        #print "trailing_y_dict found: ",trailing_y_dict      
    
        while trailing_Y_dict.keys() != []:
            max_Y_exp = -100000000
            for key in trailing_Y_dict.keys():
                if key > max_Y_exp:
                    max_Y_exp = key
            #print "right_divide_element",element
            #print "trailing_y_dict",trailing_y_dict
            P.sort()
            P.reverse()
            leading_term = P[0]
            div = trailing_Y_dict.pop(max_Y_exp)*Yi^leading_term[0]
            #print "div",div
            output = output + div
            for i in range(1,P.__len__()):
                term = P[i]
                trailing_Y_dict[max_Y_exp-leading_term[0]+term[0]] = trailing_Y_dict[max_Y_exp-leading_term[0]+term[0]] - term[1]*div*Y^term[0]
                if trailing_Y_dict[max_Y_exp-leading_term[0]+term[0]] == self.F(0):
                    trailing_Y_dict.pop(max_Y_exp-leading_term[0]+term[0])
    
            P.reverse()
          
        return output
        
    def Y_split(self, mon):
        # splits the monomial mon into its maximal non-negative Y prefix and the
        # corresponding suffix
        prefix = self.M(1)
        suffix = self.M(1)
        switch = true
        for i,j in mon._element_list:
            if i != 2 and switch:
                prefix._element_list.append([i,j])
            else:
                if j >= 0 and switch:
                    prefix._element_list.append([i,j])
                else:
                    switch = false
                    suffix._element_list.append([i,j])
                    
        return [prefix,suffix]
    