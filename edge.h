/*
 *  edge.h
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 4/25/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#ifndef edge_H
#define edge_H

class node; 

class edge {
	
public:
	edge(); //constructor
	~edge(); //destructor
	void set_pnode(node* n) {pnode = n;};
	void set_dnode(node* n) {pnode = n;};
	node* get_pnode() {return pnode;};
	node* get_dnode() {return dnode;};
	
private:
	node* pnode; //parent node
	node* dnode; //descendant node
	
};

#endif