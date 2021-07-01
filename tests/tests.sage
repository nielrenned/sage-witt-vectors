import os
import unittest

'''
Run this file with Sage from the tests directory.
'''

sage.repl.load.load('../witt.sage', globals())
magma.chdir('finotti-witt/')
magma.load('witt.m')

'''
We assume the magma code is correct. Not a bad assumption imo.
'''

class TestIntConversion(unittest.TestCase):
    def test_prime_char_conversion(self):
        for p in Primes()[:10]:
            for prec in range(1, 20):
                n = randint(0, p^prec - 1)
                with self.subTest(p=p, prec=prec, n=n):
                    W = WittRing(GF(p), prec=prec, op_method='none')
                    x = W(n)
                    y = magma.function_call('IntToWitt', [n, p, prec-1])
                    self.assertEqual(x.vec, tuple(y))
    
    @unittest.skip("We haven't implemented int conversion for non-prime characteristic bases.")
    def test_nonprime_char_conversion(self):
        pass

class TestUnaryOperations(unittest.TestCase):
    @unittest.skip("We haven't implemented negation yet.")
    def test_negation(self):
        for p in Primes()[:10]:
            for prec in range(1, 20):
                n = randint(0, p^prec - 1)
                with self.subTest(p=p, prec=prec, n=n):
                    W = WittRing(GF(p), prec=prec, op_method='none')
                    x = -W(n)
                    y = magma(f'WittNeg(IntToWitt({n1}, {p}, {prec-1}) : choice:=3)')
                    self.assertEqual(x.vec, tuple(y))
    
    @unittest.skip("We haven't implemented inversion yet.")
    def test_inversion(self):
        for p in Primes()[:10]:
            for prec in range(1, 20):
                n = randint(0, p^prec - 1)
                while n % p == 0:
                    n = randint(0, p^prec - 1)
                with self.subTest(p=p, prec=prec, n=n):
                    W = WittRing(GF(p), prec=prec, op_method='none')
                    x = 1 / W(n)
                    y = magma('WittInv(IntToWitt({n1}, {p}, {prec-1}) : choice:=3)')
                    self.assertEqual(x.vec, tuple(y))
    
    @unittest.skip("We haven't implemented inversion yet.")
    def test_inversion_of_non_unit(self):
        for p in Primes()[:10]:
            for prec in range(1, 20):
                n = randint(0, p^prec - 1)
                with self.subTest(p=p, prec=prec, n=n):
                    W = WittRing(GF(p), prec=prec, op_method='none')
                    t = list(W(n).vec)
                    t[0] = GF(p)(0)
                    with self.assertRaises(ZeroDivisionError):
                        _ = 1 / W(t)
        

# These tests should really pick a random field, and then random elements from that field
# But for now, this will do.
class TestArithmeticOperations(unittest.TestCase):
    def setUp(self):
        self.primes = list(Primes()[:6])
        self.precisions = {p : list(range(1,5)) for p in self.primes}
    
    def test_addition(self):
        for p in self.primes:
            for prec in self.precisions[p]:
                n1 = randint(0, p^prec - 1)
                n2 = randint(0, p^prec - 1)
                with self.subTest(p=p, prec=prec, n1=n1, n2=n2):
                    W = WittRing(GF(p), prec=prec, op_method='finotti_fly')
                    x = W(n1) + W(n2)
                    _ = magma.eval(f'v1 := IntToWitt({n1}, {p}, {prec-1});')
                    _ = magma.eval(f'v2 := IntToWitt({n2}, {p}, {prec-1});')
                    y = magma('WittSum(v1, v2 : choice:=3)')
                    self.assertEqual(x.vec, tuple(y))
    
    def test_subtraction(self):
        for p in self.primes:
            for prec in self.precisions[p]:
                n1 = randint(0, p^prec - 1)
                n2 = randint(0, p^prec - 1)
                with self.subTest(p=p, prec=prec, n1=n1, n2=n2):
                    W = WittRing(GF(p), prec=prec, op_method='finotti_fly')
                    x = W(n1) - W(n2)
                    _ = magma.eval(f'v1 := IntToWitt({n1}, {p}, {prec-1});')
                    _ = magma.eval(f'v2 := IntToWitt({n2}, {p}, {prec-1});')
                    y = magma('WittDiff(v1, v2 : choice:=3)')
                    self.assertEqual(x.vec, tuple(y))
    
    def test_multiplication(self):
        for p in self.primes:
            for prec in self.precisions[p]:
                n1 = randint(0, p^prec - 1)
                n2 = randint(0, p^prec - 1)
                with self.subTest(p=p, prec=prec, n1=n1, n2=n2):
                    W = WittRing(GF(p), prec=prec, op_method='finotti_fly')
                    x = W(n1) * W(n2)
                    _ = magma.eval(f'v1 := IntToWitt({n1}, {p}, {prec-1});')
                    _ = magma.eval(f'v2 := IntToWitt({n2}, {p}, {prec-1});')
                    y = magma('WittProd(v1, v2 : choice:=3)')
                    self.assertEqual(x.vec, tuple(y))
    
    @unittest.skip("We haven't implemented division yet.")
    def test_division(self):
        for p in self.primes:
            for prec in self.precisions[p]:
                n1 = randint(0, p^prec - 1)
                n2 = randint(1, p^prec - 1)
                while n2 % p == 0:
                    n2 = randint(1, p^prec - 1) 
                with self.subTest(p=p, prec=prec, n1=n1, n2=n2):
                    W = WittRing(GF(p), prec=prec, op_method='finotti_fly')
                    x = W(n1) / W(n2)
                    _ = magma.eval(f'v1 := IntToWitt({n1}, {p}, {prec-1});')
                    _ = magma.eval(f'v2 := IntToWitt({n2}, {p}, {prec-1});')
                    y = magma('WittDiv(v1, v2 : choice:=3)')
                    self.assertEqual(x.vec, tuple(y))

if __name__ == '__main__':
    unittest.main()