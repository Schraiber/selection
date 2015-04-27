/*
 *  node.h
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 4/25/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#ifndef node_H
#define node_H

#include<vector>

class edge;

class node {

public:
	node(); //constructor
	~node(void); //destructor
	void set_time(double t) {time = t;}; 
	void set_background(bool s) {sel_background = s;};
	double get_time() {return t;};
	bool get_sel_background() {return sel_background};
	virtual std::vector<edge*> get_descs();
	virtual std::vector<edge*> get_parents(); 
	
protected:
	double time;
	bool sel_background; //T if in selected background, F otherwise
};

//coalescent node
class cNode: public node {
	
public:
	void set_d1(edge* e) {d1 = e;};
	void set_d2(edge* e) {d2 = e;};
	void set_p(edge* e) {p = e;};
	edge* get_d1() {return d1;};
	edge* get_d2() {return d2;};
	edge* get_p() {return p;};
	
	std::vector<edge*> get_descs();
	std::vector<edge*> get_parents(); 
	
private:
	edge* d1; //one descendant
	edge* d2; //other descendant
	edge* p; //parent
	
};

//recombination node
class rNode: public node {

public:
	void set_p1(edge* e) {d1 = e;};
	void set_p2(edge* e) {d2 = e;};
	void set_d(edge* e) {d = e;};
	edge* get_d1() {return d1;};
	edge* get_d2() {return d2;};
	edge* get_d() {return d;};
	
	std::vector<edge*> get_descs();
	std::vector<edge*> get_parents();
	
private:
	edge* p1; //one parent
	edge* p2; //other parent
	edge* d; //descendant
};

#endif