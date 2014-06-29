#include <config.h>

#ifndef XML_H
#define XML_H

#include "utilities.h"

#define T_MAX_XML_TAG 64

typedef struct __XML_node {

  struct __XML_attr *attr;   // Pointer to the first element of a list of attributes
  int n_attr;                // Number of attributes
  struct __XML_node *next;   // Next sibling
  struct __XML_node *prev;   // Previous sibling
  struct __XML_node *parent; // Parent of this node
  struct __XML_node *child;  // Child of this node
  char *id;
  char *name;
  char *value;
}xml_node;


typedef struct __XML_attr {
  char *name;
  char *value;
  struct __XML_attr *next; // Next attribute
  struct __XML_attr *prev; // Previous attribute
}xml_attr;



void XML_Load_File(FILE *fp, xml_node *top);
int XML_Add_Character(int c, char  **bufptr, char **buffer, int *bufsize);
int XML_Parse_Element(FILE *fp, xml_node *n);
xml_node *XML_New_Node(xml_node *parent, char *name);
int XML_Set_Attribute(xml_node *n, char *attr_name, char *attr_value);
xml_attr *XML_Make_Attribute(xml_attr *prev, char *attr_name, char *attr_value);
void XML_Make_Node_Id(xml_node *n, char *id);
int XML_Set_Node_Id(xml_node *n, char *id);
int XML_Set_Node_Value(xml_node *n, char *val);
void XML_Make_Node_Value(xml_node *n, char *val);

#endif
