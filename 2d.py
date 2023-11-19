"""
Compute exponential of x iteratively
For every iteration, check among the two possible ways of doing the product and compute the best option

"""
from math import exp
import math
import matplotlib.pyplot as plt

def is_product_bounded(v,w):
	""" Returns True 
	if the product of v and w is bounded
	and False otherwise
	"""

	try:
		product = v * w
		if math.isinf(v) or v==0:
			return False
		return True
	except OverflowError:
		return False


def next_product(v,w):
	"""Returns True,v*w if the product is bounded
	and False,0 otherwise
	"""

	next = 0
	bounded = True

	if is_product_bounded(v,w):
		return True, v*w
	else:
		return False, 0


def exp_hat(x):
	""" Computes approximation of e^x iteratively
	computing the most appropriate term distribution
	for every sum term
	"""

	# Trivial case
	if x==0: return 1

	# Start with function at 0
	e_x = 1

	# Start with two options available
	opt_1 = True
	opt_2 = True
	prev_exponent = 1
	prev_factorial = 1
	bounded_x_n = True
	bounded_n_f = True
	n = 1

	# Iterative loop until no options left
	while opt_1 or opt_2:
		

		
		# Update if opt_1 and opt_2 are still available
		
		if not bounded_n_f or not bounded_x_n:
			opt_1 = False
		else:
			option1_term_1 = (prev_exponent / prev_factorial)
			option1_term_2 = (x/n)
			opt_1 = is_product_bounded(option1_term_1,option1_term_2)

		if not bounded_n_f or not bounded_x_n:
			opt_2 = False
		else:
			option2_term_1 = (prev_exponent / n)
			option2_term_2 = (x / prev_factorial)
			opt_2 = is_product_bounded(option2_term_1,option2_term_2)
		
		#opt_1 = (prev_exponent / prev_factorial) * (x/n)
		#opt_2 = (prev_exponent / n) * (x/prev_factorial)

		print(opt_1,opt_2)

		if opt_1:
			e_x += option1_term_1 * option1_term_2
		elif opt_2:
			e_x += option2_term_1 * option2_term_2
		else:
			break

		# Update terms
		#prev_exponent *= x
		#prev_factorial *= n

		# Update 4 terms if possible
		bounded_x_n, prev_exponent  = next_product(prev_exponent,x)
		bounded_n_f, prev_factorial = next_product(prev_factorial,n)
		n += 1
		
		#print(n,prev_exponent,prev_factorial)

	return e_x


print(exp_hat(10))
print(exp(10))
#print(exp_hat(0))

# Testing
test_values = [1,5,10,15,20]#[1,5,10,15,20]
e_s         = []
e_hats      = []
errors       = []

for x in test_values:
    e,e_hat = exp(x),exp_hat(x)
    error = e - e_hat
    
    e_s.append(e)
    e_hats.append(e_hat)
    errors.append(error)


# Print the values
print('x',test_values)
print('exp(x)',e_s)
print('exp_hat_2(x)',e_hats)
#print('exp(x)-exp_hat_2(x)',errors)


# Plot the functions
plt.figure(figsize=(8,2),dpi=400)
plt.title('Comparisson of exp(x) and our exp_hat_2(x) function')
plt.semilogy(test_values,e_s,label='exp(x)')
plt.semilogy(test_values,e_hats,label='exp_hat(x)')
plt.xlim([-20,-1])
plt.xlabel ('x')
plt.ylabel ('y')
plt.ylim([e_s[0],e_s[len(test_values)-1]])
plt.legend()
plt.tight_layout()
plt.show()

# Plot the errors
plt.figure(figsize=(8,2),dpi=400)
plt.title('Comparisson of exp(x) and exp_hat_2(x) errors')
plt.plot(test_values,errors,label='Error',lw=2,c='green')
plt.xlim([-20,-1])
plt.xlabel ('x')
plt.ylabel ('exp(x)-exp_hat_2(x)')
plt.legend()
plt.tight_layout()
plt.show()
