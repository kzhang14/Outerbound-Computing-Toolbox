# computeOuterbound

This program takes user-defined entropies of communication system or data storage system, and outputs the outerbound of system as information inequalities.

Instructions:

1. Convert the problem descriptions to mathematical equalities/inequalities, usually in the format of joint entropy terms or mutual information terms, and put them into the 'input.txt' file, the file name can be arbitrary.

2. The 'input.txt' file has certain format:

* It has 3 sections: Variables, Minimize (or Maximize), Subject To, which together defines the minimum possible structure of a constrained **Linear Programming** problem.

* 'Variable' section:
** The 'Variable' needs to be only-first-letter capitalized, the ending colon can be omitted. 
** Starting a new line, define the variable names that's going to be used to describe the system. Variable name can include letters, numbers and the following symbols: !@#$%^&*~.;`", but can not include {}{}<>()|=-+/\,:, nor the  following letters: H, I (which are reserved to denote joint entropy and mutual information).

* 'Minimize' section:
** Can use 'Maximize'
** Starting a new line, describe the objective function needs to be optimized, as a function of joint entropy or mutual information entities.

* 'Subject To' section:
** Each new line defines a equality or inequality function of joint entropy or mutual information entities.


3. The symmetry of the system can be defined in another file, 'perm.txt', usually in the format of a symmetry group. Each line of the file denotes a permutation in the group (mathematical concept), in the cycle format.  


4. Call the function ```formatConvert.py``` with input file as argument 
```
>> python formatConvert.py input.txt
```

The function will take the input.txt and convert to formal **MPS** or **LP** file format which can be further used by other optimization softwares, such as **IBM CPLEX** and **Gurobi**. 




