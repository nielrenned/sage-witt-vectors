# Raise x^n power in characteristic p.
# If x is not an element of a ring of characteristic p, throw an error.
# If x is an element of GF(p^k), this is already fast.
# However, if x is a polynomial, we need to do it better (wtf Sage).
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_element import Polynomial

def _fast_char_p_power(x, n):
    if n not in ZZ:
        raise ValueError(f'Exponent {n} is not an integer')
    if n == 0 or x == 1:
        return x.parent().one()
    if x.parent().characteristic() not in Primes():
        raise ValueError(f'{x} is not in a ring of prime characteristic')
    
    is_Polynomial  = isinstance(x, Polynomial)
    is_MPolynomial = isinstance(x, MPolynomial)
    
    if not (is_Polynomial or is_MPolynomial):
        return x^n
    if (is_Polynomial and x.is_gen()) or (is_MPolynomial and x.is_generator()):
        return x^n
    if n < 0:
        x = x^-1 # This may throw an error.
        n = -n
    
    P = x.parent()
    p = P.characteristic()
    base_p_digits = ZZ(n).digits(base=p)
    
    x_to_the_n = 1
    
    for p_exp, digit in enumerate(base_p_digits):
        if digit == 0:
            continue
        inner_term = x^digit
        term_dict = {}
        for e_int_or_tuple, c in inner_term.dict().items():
            power = p^p_exp
            new_c = _fast_char_p_power(c, power)
            #new_e_tuple = e_tuple.emul(power)
            new_e_tuple = None
            if is_Polynomial: # Then the dict keys are ints
                new_e_tuple = e_int_or_tuple * power
            elif is_MPolynomial: # Then the dict keys are ETuples
                new_e_tuple = e_int_or_tuple.emul(power)
            term_dict[new_e_tuple] = new_c
        term = P(term_dict)
        x_to_the_n *= term
    
    return x_to_the_n

_fcppow = _fast_char_p_power

class WittVector(CommutativeRingElement):
    def __init__(self, parent, vec=None):
        self.prec = parent.precision()
        B = parent.base()
        if vec is not None:
            if len(vec) != self.prec:
                raise ValueError(f'{vec} is not the correct length. Expected length to be {self.prec}.')
            self.vec = tuple(B(x) for x in vec)
        else:
            self.vec = (B(0) for i in range(self.prec))
        CommutativeRingElement.__init__(self, parent)
    
    def _richcmp_(self, other, op):
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            return self.vec == other.vec
        elif op == op_NE:
            return self.vec != other.vec
        else:
            return NotImplemented
    
    def _repr_(self):
        return '(' + ', '.join(map(str, self.vec)) + ', ...)'
    
    def _add_(self, other):
        P = self.parent()
        C = self.__class__
        # As a slight optimization, we'll check for zero ahead of time.
        # This also has the benefit of allowing us to create polynomials, even
        # if we choose 'none' as the op_method.
        if other == P.zero():
            return self
        elif self == P.zero():
            return other
        
        op_method = P._op_method()
        if op_method == 'polys':
            s = P.sum_polys
            sum_vec = tuple(s[i](*(self.vec + other.vec)) for i in range(self.prec)) # note here this is tuple addition, i.e. concatenation
            return C(P, vec=sum_vec)
        elif op_method.startswith('finotti'):
            x = self.vec
            y = other.vec
            prec = P.precision()
            p = P.prime
            bin_table = P.binomial_table
            G = []
            for n in range(0, prec):
                G_n = [x[n], y[n]]
                for i in range(0, n):
                    G_n.append(P.eta_bar(G[i], n-i))
                G.append(G_n)
            sum_vec = tuple(sum(G[i]) for i in range(prec))
            return C(P, vec=sum_vec)
        else:
            return NotImplemented
    
    def _mul_(self, other):
        P = self.parent()
        C = self.__class__
        # As a slight optimization, we'll check for zero or one ahead of time.
        # This also has the benefit of allowing us to create polynomials, even
        # if we choose 'none' as the op_method
        if other == P.zero():
            return P.zero()
        elif self == P.zero():
            return P.zero()
        elif other == P.one():
            return self
        elif self == P.one():
            return other
        
        op_method = P._op_method()
        if op_method == 'polys':
            # note here this is tuple addition, i.e. concatenation
            p = P.prod_polys
            prod_vec = tuple(p[i](*(self.vec + other.vec)) for i in range(self.prec))
            return C(P, vec=prod_vec)
        elif op_method.startswith('finotti'):
            x = self.vec
            y = other.vec
            G = [[x[0] * y[0]]]
            prec = P.precision()
            p = P.prime
            bin_table = P.binomial_table
            for n in range(1, prec):
                #G_n = [x[0]^(p^n) * y[n], y[0]^(p^n) * x[n]]
                #G_n.extend(x[i]^(p^i) * y[n-i]^(p^(n-i)) for i in range(1,n))
                G_n = [_fcppow(x[0], p^n) * y[n], _fcppow(y[0], p^n) * x[n]]
                G_n.extend(_fcppow(x[i], p^i) * _fcppow(y[n-i], p^(n-i)) for i in range(1, n))
                for i in range(0, n):
                    G_n.append(P.eta_bar(G[i], n-i))
                G.append(G_n)
            prod_vec = tuple(sum(G[i]) for i in range(prec))
            return C(P, vec=prod_vec)
        else:
            return NotImplemented
    
    def __invert__(self):
        if self.vec[0].is_zero():
            raise ZeroDivisionError(f"Inverse of {self} does not exist.")
        P = self.parent()
        if self.prec == 1:
            return P((self.vec[0]^-1, ))
        var_names = [f'Y{i}' for i in range(1, self.prec)]
        poly_ring = PolynomialRing(P.base(), var_names)
        gens = poly_ring.gens()
        inv_vec = list( (self.vec[0]^-1, ) + gens )
        W = WittRing(poly_ring, prec=P.prec, op_method=P.op_method)
        prod_vec = (W(self.vec)*W(inv_vec)).vec
        for i in range(1, self.prec):
            poly = prod_vec[i](inv_vec[1:])
            inv_vec[i] = -poly.constant_coefficient() / poly.monomial_coefficient(gens[i-1])
        return P(inv_vec)
    
    def __neg__(self):
        P = self.parent()
        if P.prime == 2:
            all_ones = P(tuple(1 for _ in range(self.prec)))
            return all_ones*self
        neg_vec = tuple(-self.vec[i] for i in range(self.prec))
        return P(neg_vec)
    
    def xi(self, n):
        P = self.parent()
        p = P.prime
        return P.teichmuller_lift(self.vec[n].nth_root(p^n))
    
    def proj(self):
        return self.vec[0]
    
    def proj_xi(self, n):
        p = self.parent().prime
        return self.vec[n].nth_root(p^n)
    
    def frobenius(self, n=1):
        P = self.parent()
        C = self.__class__
        p = P.prime
        return C(P, vec=tuple(_fcppow(x, p^n) for x in self.vec))

# We _should_ convert to an extension of Zp if the base is a finite field of prime characteristic.
# But the point is to compare, so this is better for now
class WittRing(CommutativeRing, UniqueRepresentation):
    Element = WittVector
    
    # The options for operation_method: 'none', 'polys', 'finotti_fly', 'finotti_pregen'
    # The 'polys' option is _incredibly_ slow.
    # Note that 'prec' is the length of the internal representation. (At a cursory
    # glance, this appears to agree with the 'n' from the paper. It probably doesn't matter.)
    # If prime is not supplied, we attempt to use the characteristic of the base. If this is not prime,
    # we throw an error.
    def __init__(self, base, prec=4, prime=None, category=None, op_method='none'):
        self.prec = prec
        
        # In theory, we can create W_n(R) for any commutative ring R, but we don't know how to
        # convert integers to W_n(R) without just adding 1 to itself over and over.
        # This makes it very hard to work with, so for now, we'll just error until we get a better idea
        # for what to do.
        if base not in CommutativeRings():
            raise ValueError(f'Cannot create Ring of Witt Vectors over {base}: {base} is not a commutative ring.')
        
        if prime is None:
            char = base.characteristic()
            if char not in Primes():
                raise ValueError(f'Cannot create Ring of Witt Vectors over {base}: {base} has non-prime characteristic, and no prime was supplied.')
            else:
                self.prime = char
        else:
            self.prime = prime
        
        if op_method not in ['none', 'polys', 'finotti_fly', 'finotti_pregen']:
            raise ValueError(f'"{op_method}" is not a valid operation method for Witt Rings. It must be one of ["none", "polys", "finotti_fly", "finotti_pregen"]')
        
        if self.prime != base.characteristic() and op_method.startswith('finotti'):
            raise ValueError(f'Currently, the finotti methods for sum and product don\'t work if char(base) != p.')
        
        '''
        if base not in CommutativeRings() or base.characteristic() not in Primes():
            raise ValueError(f'Can only create Witt Rings over commutative rings of prime characteristic.')
        
        self.prime = base.characteristic()
        '''
        
        self.op_method = op_method
        self.sum_polys = None
        self.prod_polys = None
        self.eta_bars = None
        self.eta_bar_cache = {}
        self.greenberg_transform_Ds_cache = {}
        
        if op_method == 'polys':
            self.generate_sum_and_product_polynomials(base)
        elif op_method.startswith('finotti'):
            self.generate_binomial_table()
            if op_method == 'finotti_pregen':
                self.generate_eta_bars(base)
            else:
                pass # don't need to do anything else
        elif op_method == 'none':
            pass # do nothing
        
        if category is not None:
            if base.is_integral_domain():
                category = IntegralDomains()
            else:
                category = CommutativeRings()
        CommutativeRing.__init__(self, base, category=category)
    
    def _repr_(self):
        return f"Ring of Witt Vectors of length {self.prec} over {self.base()}"
    
    def characteristic(self):
        return 0
    
    def precision(self):
        return self.prec
    
    def _coerce_map_from_(self, S):
        if S is ZZ or S is GF(self.prime):
            return True
        else:
            return False
    
    def _element_constructor_(self, x):
        if x in ZZ:
            return self.element_class(self, self.int_to_vector(x))
        elif isinstance(x, tuple) or isinstance(x, list):
            return self.element_class(self, x)
        else:
            return NotImplemented
    
    def random_element(self):
        return self.element_class(self, tuple(self.base().random_element() for _ in range(self.prec)))
    
    def generate_sum_and_product_polynomials(self, base):
        p = self.prime
        prec = self.prec
        x_var_names = ['X{}'.format(i) for i in range(prec)]
        y_var_names = ['Y{}'.format(i) for i in range(prec)]
        var_names = x_var_names + y_var_names
        # Okay, wtf. Sage, by default relies on Singular for Multivariate Polynomial Rings, 
        # but Singular uses only SIXTEEN bits (unsigned) to store its exponents. So if we 
        # want exponents larger than 2^16 - 1, we have to use the generic implementation.
        # However, after some experimentation, it seems like the generic implementation
        # is faster? Idk, I'm confused.
        # After trying to compute S_4 for p=5, it looks like generic is faster for 
        # very small polys, and MUCH slower for large polys. So we'll default to singular
        # unless we can't use it.
        # 
        # Remark: Since when is SIXTEEN bits sufficient for anyone???
        #
        if p^(prec-1) >= 2^16:
            implementation = 'generic'
        else:
            implementation = 'singular'
        R = PolynomialRing(ZZ, var_names, sparse=True, implementation=implementation)
        x_y_vars = R.gens()
        x_vars = x_y_vars[:prec]
        y_vars = x_y_vars[prec:]
        
        self.sum_polys = {}
        for n in range(prec):
            s_n = x_vars[n] + y_vars[n]
            for i in range(n):
                s_n += (x_vars[i]^(p^(n-i)) + y_vars[i]^(p^(n-i)) - self.sum_polys[i]^(p^(n-i))) / p^(n-i)
            self.sum_polys[n] = R(s_n) # Dividing by p^k converts it to a polynomial over Q, so convert back to Z
        
        self.prod_polys = {0: x_vars[0] * y_vars[0]}
        for n in range(1, prec):
            x_poly = sum([p^i * x_vars[i]^(p^(n-i)) for i in range(n+1)])
            y_poly = sum([p^i * y_vars[i]^(p^(n-i)) for i in range(n+1)])
            p_poly = sum([p^i * self.prod_polys[i]^(p^(n-i)) for i in range(n)])
            p_n = (x_poly*y_poly - p_poly) / p^n
            self.prod_polys[n] = R(p_n) # Dividing by p^k converts it to a polynomial over Q, so convert back
        
        # We have to use generic here, because Singular doesn't support Polynomial Rings over Polynomial Rings
        # e.g. PolynomialRing(GF(5)['x'], ['X', 'Y'], implementation='singular') will fail
        S = PolynomialRing(base, x_y_vars, sparse=True, implementation='generic')
        for n in range(prec):
            self.sum_polys[n] = S(self.sum_polys[n])
            self.prod_polys[n] = S(self.prod_polys[n])
    
    def generate_binomial_table(self):
        import numpy as np
        p = self.prime
        R = Zp(p, prec=self.prec+1, type='fixed-mod')
        v_p = ZZ.valuation(p)
        table = [[0]]
        for k in range(1, self.prec+1):
            row = np.empty(p^k, dtype=int)
            row[0] = 0
            prev_bin = 1
            for i in range(1, p^k // 2 + 1):
                val = v_p(i)
                # Instead of calling binomial each time, we compute the coefficients
                # recursively. This is MUCH faster.
                next_bin = prev_bin * (p^k - (i-1)) // i 
                prev_bin = next_bin
                series = R(-next_bin // p^(k-val))
                for _ in range(val):
                    temp = series % p
                    series = (series - R.teichmuller(temp)) // p
                row[i] = ZZ(series % p)
                row[p^k - i] = row[i] # binomial coefficients are symmetric
            table.append(row)
        self.binomial_table = table
    
    def generate_eta_bars(self, base):
        # Here we're going to generate the two-variable eta-bars ahead of time, that
        # way, we don't have to recompute over and over.
        R.<x,y> = base['x', 'y'] # TODO: This should maybe be F_p right? That would make precomputing much faster.
        self.eta_bars = {}
        self.eta_bars[0] = x+y # Technically, we don't need this. It's included for completeness.
        for k in range(1, self.prec+1):
            self.eta_bars[k] = self.eta_bar((x, y), k)
    
    def eta_bar(self, vec, eta_index):
        #res_str = f'eta_{eta_index}({tuple(vec)}) ='
        vec = tuple(x for x in vec if x != 0) # strip zeroes
        
        # special cases
        ret_val = None
        if len(vec) <= 1:
            #ret_val = 0
            return 0
        if ret_val is None and eta_index == 0:
            ret_val = sum(vec)
        
        # Check cache.
        # We use the frozenset of vec after stripping zeroes, since order doesn't
        # matter and we can remove zeroes (by Remarks 5.2).
        cache_key = (frozenset(vec), eta_index)
        if cache_key in self.eta_bar_cache:
            ret_val = self.eta_bar_cache[cache_key]
            #res_str = '*' + res_str
            #return self.eta_bar_cache[cache_key]
        
        # compute bin_table if necessary
        #if bin_table is None:
        #    bin_table = GT_Bin_Table(prime, eta_index)
        
        # renaming to match notation in paper
        k = eta_index
        p = self.prime
        # if vec = (x,y), we know what to do: Theorem 8.6
        if ret_val is None and len(vec) == 2:
            # Here we have to check if we've pre-computed already
            if self.op_method == 'finotti_pregen' and (k in self.eta_bars):
                return self.eta_bars[k](*vec)
            else:
                x,y = vec
                scriptN = [[None] for _ in range(k+1)] # each list starts with None, so that indexing matches paper
                # calculate first N_t scriptN's
                for t in range(1, k+1):
                    for i in range(1, p^t):
                        scriptN[t].append(self.binomial_table[t][i] * _fcppow(x, i) * _fcppow(y, p^t - i))
                indexN = [p^i - 1 for i in range(k+1)]
                for t in range(2, k+1):
                    for l in range(1,t):
                        # append scriptN_{t, N_t+l}
                        next_scriptN = self.eta_bar(scriptN[t-l][1:indexN[t-l]+t-l], l)
                        scriptN[t].append(next_scriptN)
                #ret_val = sum(scriptN[k][1:])
                return sum(scriptN[k][1:])
        
        # if vec is longer, we split and recurse: Proposition 5.4
        # This is where we need to using multiprocessing.
        if ret_val is None and len(vec) > 2:
            m = len(vec) // 2
            v_1 = vec[:m]
            v_2 = vec[m:]
            s_1 = sum(v_1)
            s_2 = sum(v_2)
            total = 0
            scriptM = [[] for _ in range(k+1)]
            for t in range(1, k+1):
                scriptM[t].append(self.eta_bar(v_1,        t))
                scriptM[t].append(self.eta_bar(v_2,        t))
                scriptM[t].append(self.eta_bar((s_1, s_2), t))
            for t in range(2, k+1):
                for s in range(1, t):
                    result = self.eta_bar(scriptM[t-s], s)
                    scriptM[t].append(result)
            #ret_val = sum(scriptM[k])
            return sum(scriptM[k])
        
        # Cache and return
        self.eta_bar_cache[cache_key] = ret_val
        #if verbose:
        #    print(res_str, ret_val)
        return ret_val
    
    def _op_method(self):
        return self.op_method
    
    # Currently, this only works if the base has characteristic p.
    # An open question (I think) is: Can we find an efficient method to embed Z in the WittRing with _any_ base?
    def int_to_vector(self, n):
        p = self.prime
        R = Zp(p, prec=self.prec+1, type='fixed-mod')
        F = GF(p)
        B = self.base()
        series = R(n)
        witt_vector = []
        for _ in range(self.prec):
            # Probably slightly faster to do "series % p," but this way, temp is in F_p
            temp = F(series)
            witt_vector.append(B(temp)) # make sure elements of vector are in base
            series = (series - R.teichmuller(temp)) // p
        return witt_vector
    
    def teichmuller_lift(self, x):
        if x not in self.base():
            raise Exception(f'{x} not in {self.base()}')
        else:
            return self.element_class(self, (x,) + tuple(0 for _ in range(self.prec-1)))
    
    def greenberg_transform(self, f):
        if f.parent().base() is not self or len(f.parent().gens()) != 2:
            raise ValueError(f'f must be a polynomial in x,y over {self}')
        p = self.prime
        G = []
        for n in range(self.prec):
            f_sigma = f.map_coefficients(lambda c: c.frobenius(n))
            G_n = []
            for (x_e, y_e), c in zip(f_sigma.exponents(), f_sigma.coefficients()):
                sigma_c = c.frobenius(n)
                for r in range(0, n+1):
                    for i in range(0, r+1):
                        D = self.greenberg_transform_Ds(i, r-i, n)
                        S = D[r][0].parent()
                        x0, y0 = S.gens()[0], S.gens()[self.prec]
                        if x_e - i < 0 or y_e - (r-i) < 0:
                            continue
                        deriv_c = sigma_c * binomial(x_e, i) * binomial(y_e, r-i)
                        for j in range(r, n+1):
                            xi_c = deriv_c.proj_xi(n-j)
                            if xi_c == 0:
                                continue
                            for k in range(r, j+1):
                                if k not in D or D[k][j-k] == 0:
                                    continue
                                G_n.append(xi_c * _fcppow(x0, ((x_e-i)*p^n)) * _fcppow(y0, ((y_e-(r-i))*p^n)) * D[k][j-k])
            for i in range(n):
                G_n.append(self.eta_bar(G[i], n-i))
            G.append(G_n)
        return tuple(map(sum, G))
    
    # Return D^{i,j}_{k,n,l} for each k in {i+j, ..., n(i+j)} and l in {0, ..., n-1} 
    # Format is a dict of lists indexed by [k][l]
    def greenberg_transform_Ds(self, i, j, n):
        if (i, j, n) in self.greenberg_transform_Ds_cache:
            return self.greenberg_transform_Ds_cache[(i,j,n)]
        p = self.prime
        D = {}
        x_var_names = ['x{}'.format(i) for i in range(self.prec)]
        y_var_names = ['y{}'.format(i) for i in range(self.prec)]
        var_names = x_var_names + y_var_names
        R = PolynomialRing(Integers(p^(n+1)), var_names, sparse=True)
        x_y_vars = R.gens()
        x_vars = x_y_vars[:self.prec]
        y_vars = x_y_vars[self.prec:]
        S.<t> = R['t']
        x_expr = S(0)
        y_expr = S(0)
        for e in range(1, n+1):
            x_expr += t^(e) * x_vars[e]^(p^(n-e))
            y_expr += t^(e) * y_vars[e]^(p^(n-e))
        f = x_expr^i * y_expr^j
        
        W = WittRing(GF(p), prec=n+1, op_method='none')
        for k in range(i+j, n*(i+j) + 1):
            coeff = f.monomial_coefficient(t^k)
            D_k = [0 for _ in range(n+1)]
            for mon in coeff.monomials():
                witt_vec = W(coeff.monomial_coefficient(mon))
                for l in range(n+1):
                    D_k[l] += witt_vec.proj_xi(l)*mon
            D[k] = D_k
        self.greenberg_transform_Ds_cache[(i,j,n)] = D
        return D