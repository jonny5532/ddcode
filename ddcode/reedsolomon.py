from itertools import izip_longest
import math

"""Basic Reed-Solomon error correction operations for PDF417 barcode data."""

def invert(a, p):
	"""Invert an integer a in GF(p)"""
	u,v = a,p
	x1, x2 = 1, 0
	while u!=1:
		q = v/u
		r, x = v-q*u , x2-q*x1
		v,u = u, r
		x2,x1 = x1,x
	return x1%p

class GFPoly:
	"""Represent and perform basic operations on a polynomial in GF929.
	
	Basic usage:
	p = GFPoly([1,2,3]) #1 is x^0 term, 2 is x^1 term, etc.
	
	"""
	
	def __init__(self, poly=[0], modulo=929):
		self.poly = poly #list of coefficients starting at x^0
		while len(self.poly) and self.poly[-1]==0: self.poly = self.poly[:-1]
		self.modulo = modulo
		
	def degree(self): #returns highest power in poly
		return len(self.poly)-1
	
	def invert(self):
		assert self.degree()==0 #only works on constants atm
		return GFPoly([invert(self.poly[0], self.modulo)], self.modulo)
	
	def mod(self, g):
		remainder = [0]*(max(self.degree(), g.degree())+1)
		for i,t in list(enumerate(self.poly))[::-1]:
			deg = i-g.degree()
			c = ((t/g.poly[-1]) - remainder[deg+g.degree()])%self.modulo
			remainder[deg+g.degree()] = 0
			for j,t in list(enumerate(g.poly[:-1]))[::-1]:
				remainder[j+deg] = (remainder[j+deg] + (c*t))%self.modulo
			if deg==0: break
		return GFPoly([(-i)%929 for i in remainder], modulo=self.modulo)

	def truncate(self, pwr): #truncate orders pwr and above
		return GFPoly(self.poly[:pwr])

	def evaluate(self, x):
		x = x%929
		xs = [p*(x**i) for i,p in enumerate(self.poly)]
		return sum(xs)%self.modulo

	def __getitem__(self, k):
		if k>=len(self.poly): return 0
		return self.poly[k]

	def __setitem__(self, k, v):
		while len(self.poly)<=k: self.poly.append(0)
		self.poly[k] = v%self.modulo

	def __repr__(self):
		ret = ""
		for i,c in list(enumerate(self.poly))[::-1]:
			if i==0:
				ret += "%d"%c
			elif c==0: continue
			else:
				eval
				ret += "%dx%s+"%(c,unichr(63498+i))
		return ret.encode("utf-8")
		
	def __sub__(self, p):
		return GFPoly([(a-b)%self.modulo for a,b in izip_longest(self.poly, p.poly, fillvalue=0)])
		
	def __add__(self, p):
		return GFPoly([(a+b)%self.modulo for a,b in izip_longest(self.poly, p.poly, fillvalue=0)])
		
	def __mul__(self, p):
		#very naiive multiplication!
		ret = GFPoly()
		for i,a in enumerate(self.poly):
			if a==0: continue
			for j,b in enumerate(p.poly):
				if b!=0: ret[i+j] = ret[i+j] + (a*b)%929
		return ret
		
	def berlekamp_massey(self, datalen, t):
		s = [self.evaluate(3**i) for i in range(1, t+1)]
		#print "syndromes are",s

		C,B,x,L,b,N = GFPoly([1]), GFPoly([1]), 1, 0, 1, 0

		while N<t:
			d = (s[N] + sum([C[i]*s[N-i] for i in range(1, L+1)]))%929
			if d==0:
				pass
			elif 2*L > N:
				t1 = GFPoly()
				t1[x] = d * invert(b, 929)
				C = C - t1*B
				x = x + 1
			else:
				T = C
				t1 = GFPoly()
				t1[x] = d * invert(b, 929)
				C = C - t1*B
				L = N + 1 - L
				B = T
				b = d
				x = 1
			N = N + 1
			

		roots = []
		for i in range(929):
			if C.evaluate(i)==0: roots.append(i)
			
		out_roots, error_positions = [], []
		
		for pos in range(datalen):
			inv = invert((3**pos), 929)
			if inv in roots:
				out_roots.append(inv)
				error_positions.append(pos)

		#print "error positions:",error_positions

		return GFPoly(s),C,out_roots,error_positions

		
	def formal_derivative(self):
		return GFPoly([(self.poly[i]*i)%929 for i,v in enumerate(self.poly)][1:])


def pdf417_rs_encode(data, poly):
	pp = [0]*(len(poly)) + data[::-1] 
	g = poly + [1]
	
	leftovers = [0]*len(pp)

	for i,t in list(enumerate(pp))[::-1]:
		deg = i-len(g)+1
		c = ((t/g[-1]) - leftovers[deg+len(g)-1])%929
		leftovers[deg+len(g)-1] = 0
		for j,t in list(enumerate(g[:-1]))[::-1]:
			leftovers[j+deg] = (leftovers[j+deg] + (c*t)%929)%929
		if deg==0: break
	s = [(x+y) for x,y in izip_longest(pp, leftovers, fillvalue=0)]
	return s[::-1]

alpha = 3 #seems to be a PDF417 feature

def pdf417_rs_decode(s, t):
	s = [(0 if i is None else i) for i in s]
	
	r = GFPoly(s[::-1])

	S, C, error_roots, error_positions = r.berlekamp_massey(len(s), t)
	sigma = (S*C).truncate(t)#.truncate(2*len(error_positions)) #4 being 2t, where t is number of errors?
	fd = C.formal_derivative()

	for er,ep in zip(error_roots, error_positions):
		r[ep] = r[ep] - ((-sigma.evaluate(er))%929)*invert(fd.evaluate(er), 929)%929

	return r.poly[::-1]


if __name__=="__main__":
	#silly example
	data = [3,2,1]
	rspoly = [522, 568, 723, 809]
	print "Encoding",data,"with generator",rspoly
	s = pdf417_rs_encode(data, rspoly)
	print "Message is",s
	
	d = len(rspoly)+1
	a = 3
	
	s[2] = None
	s[3] = 7
	
	print "Damaged message is",s
	print "Corrected message is",pdf417_rs_decode(s, 4)
