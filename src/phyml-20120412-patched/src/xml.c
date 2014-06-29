#include "xml.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Load_File(FILE *fp, xml_node *top)
{
  int c;
  char *buffer,*bufptr;
  int bufsize;
  xml_node *parent;
  xml_node *node;

  buffer = (char *)mCalloc(T_MAX_XML_TAG,sizeof(char));

  bufsize = T_MAX_XML_TAG;
  bufptr  = buffer;
  parent  = top;
  node    = NULL;

  while((c = fgetc(fp)) != EOF)
    {
      if(c == '<' && bufptr > buffer) 
	{
	  *bufptr = '\0';

	  PhyML_Printf("\n. Read value '%s' for node '%s'",buffer,node->name);
	  fflush(NULL);
	  XML_Set_Node_Value(node,buffer);
	  bufptr = buffer;
	}
	      
      if(c == '<')
	{
	  bufptr = buffer;

	  while((c = fgetc(fp)) != EOF)
	    {
	      if(isspace(c) != NO || c == '>' || (c == '/' && bufptr > buffer)) break; // End of open or close tag
	      else if(c == '<')
		{
		  Exit("\n. Bare < in element!");
		}	      
	      else if(XML_Add_Character(c,&bufptr,&buffer,&bufsize))
		{
		  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("\n");	  
		}
	    }

	  *bufptr = '\0';
	  
	  if(!strcmp(buffer,"!--")) // Get the rest of the comment
	    {
	      while((c = fgetc(fp)) != EOF)
		{
		  
		  if(c == '>' && bufptr > (buffer + 4) && bufptr[-3] != '-' &&
		     bufptr[-2] == '-' && bufptr[-1] == '-') break;
		  else if(XML_Add_Character(c,&bufptr,&buffer,&bufsize))
		    {
		      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		      Exit("\n");	  
		    }
		}
	      *bufptr = '\0';

	      if(c != '>')
		{
		  PhyML_Printf("\n. Early EOF in comment node.");
		  Exit("\n");	  
		}	      
	    }	  
	  else if(buffer[0] == '/') // Close tag
	    {
	      if(strcmp(buffer+1,parent->name))
		{
		  PhyML_Printf("\n. Opened tag with name '%s' and closed it with '%s'...",node->name,buffer+1);
		  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("\n");
		}

	      printf("\n. Closing node with name '%s'",node->name);
	      parent = parent->parent;
	      node   = parent;
	    }
	  else if(buffer[0] == '?')
	    {
	      while((c = fgetc(fp)) != EOF)
		{
		  if (c == '>' && bufptr > buffer && bufptr[-1] == '?')
		    break;
		  else if (XML_Add_Character(c, &bufptr, &buffer, &bufsize))
		    {
		      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		      Exit("\n");	  
		    }
		}

	      if(c != '>')
		{
		  PhyML_Printf("\n. An error occurred when reading the processing instruction.");
		  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("\n");	  
		}

	      *bufptr = '\0';

	    }
	  else // Open tag
	    {
	      printf("\n. Create node with name '%s'",buffer);
	      node = XML_New_Node(parent,buffer);
	      if(isspace(c) != NO) c=XML_Parse_Element(fp,node);
	      else if(c == '/')
		{
		  if((c=fgetc(fp)) != '>')
		    {
		      PhyML_Printf("\n. Expected '>' but read '%c' instead",c);
		      Exit("\n");
		    }
		  c = '/';
		}

	      if(c != '/') parent = node;

	      buffer[0] = '\0';
	    }	  
	  bufptr = buffer;
	}
      else if(isspace(c) == NO)
	{
	  if(XML_Add_Character(c,&bufptr,&buffer,&bufsize))
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("\n");
	    }
	} 
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Add_Character(int c, char  **bufptr, char **buffer, int *bufsize)
{
  char *newbuffer;

  if(*bufptr >= (*buffer + *bufsize - 4))
    {
      // Increase the size of the buffer...
      
      if (*bufsize < 1024)
	(*bufsize) *= 2;
      else
	(*bufsize) += 1024;

    if ((newbuffer = realloc(*buffer, *bufsize)) == NULL)
    {
      Free(*buffer);
      PhyML_Printf("Unable to expand string buffer to %d bytes!", *bufsize);
      Exit("\n");
    }

    *bufptr = newbuffer + (*bufptr - *buffer);
    *buffer = newbuffer;
  }

  *(*bufptr)++ = c;
  return 0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Parse_Element(FILE *fp, xml_node *n)
{
  int c;
  int quote;
  char *name, *value, *ptr;
  int namesize, valsize;

  name  = (char *)mCalloc(64,sizeof(char));
  value = (char *)mCalloc(64,sizeof(char));
  
  namesize = 64;
  valsize  = 64;
  
  while((c = fgetc(fp)) != EOF)
    {

      if(isspace(c) != NO) continue;

      if(c == '/') // End of tag
	{
	  printf("\n. Closing node '%s'.",n->name);

	  quote = fgetc(fp);
	  if(quote != '>')
	    {
	      PhyML_Printf("\n. Expected '>' after '%c' but read '%c' instead",c,quote);
	      Exit("\n");
	    }
	  break;
	}
      else if(c == '<')
	{
	  Exit("\n. Bare < in element!");	  
	}
      else if(c == '>') // End of tag
	{
	  break;
	}

      name[0] = c;
      ptr     = name + 1;

      if(c == '\"' || c == '\'') // Name is in quotes
	{
	  quote = c;

	  while((c = fgetc(fp)) != EOF)
	    {
	      if(XML_Add_Character(c,&ptr,&name,&namesize))
		{
		  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("\n");
		}
	      if(c == quote) break;
	    }	  
	}
      else // Name not in quotes
	{
	  while((c = fgetc(fp)) != EOF)
	    {
	      if(isspace(c) != NO || c == '=' || c == '/' || c == '>' || c == '?')
		break;
	      else
		{
		  if(XML_Add_Character(c,&ptr,&name,&namesize))
		    {
		      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		      Exit("\n");
		    }
		}	  
	    }
	}
      
      *ptr = '\0';
            
      while(c != EOF && isspace(c) != NO) c = fgetc(fp);

      if(c == '=') // Read the attribute value
	{
	  while((c = fgetc(fp)) != EOF && isspace(c) != NO);

	  if(c == EOF)
	    {
	      PhyML_Printf("\n. Missing value in attribute.");
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("\n");
	    }

	  if(c == '\'' || c == '\"')
	    {
	      quote = c;
	      ptr   = value;

	      while((c = fgetc(fp)) != EOF)
		{
		  if(c == quote) break;
		  else
		    {
		      if(XML_Add_Character(c,&ptr,&value,&valsize))
			{
			  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			  Exit("\n");
			}
		    }
		}
	      *ptr = '\0';
	    }
	  else
	    {
	      value[0] = c;
	      ptr      = value + 1;
	      
	      while((c = fgetc(fp)) != EOF)
		{
		  if(isspace(c) != NO || c == '=' || c == '/' || c == '>')
		    break;
		  else
		    {
		      if(XML_Add_Character(c,&ptr,&value,&valsize))
			{
			  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			  Exit("\n");
			}		      
		    }
		}	      
	    }
	}

      printf("\n. Setting attribute '%s=%s' to node '%s'",name,value,n->name);
      XML_Set_Attribute(n,name,value);

      if(c == '>') break;


    }
  Free(name);
  Free(value);

  printf("\n. Return '%c'\n",c);
  return(c);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_node *XML_New_Node(xml_node *parent, char *name)
{
  xml_node *new_node;

  new_node = (xml_node *)mCalloc(1,sizeof(xml_node));

  if(name)
    {
      new_node->name = (char *)mCalloc(strlen(name)+1,sizeof(char));
      strcpy(new_node->name,name);
    }

  new_node->parent =  parent ? parent : NULL;

  new_node->next   = NULL;
  new_node->prev   = NULL;
  new_node->child  = NULL;

  if(parent)
    { 
      xml_node *prev;

      prev = NULL;

      while(parent->child) 
	{
	  prev          = parent->child;
	  parent->child = parent->child->next;
	}
      parent->child = new_node;
      new_node->prev = prev;
      if(prev) new_node->prev->next = new_node;
    }

  new_node->attr   = NULL;
  
  return new_node;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Set_Attribute(xml_node *n, char *attr_name, char *attr_value)
{
  xml_attr *prev;
  char *s;

  prev = NULL;
  while(n->attr != NULL) 
    {
      n->attr = n->attr->next;
      prev    = n->attr;
    }

  n->attr = XML_Make_Attribute(prev,attr_name,attr_value);
  n->n_attr++;

  while(n->attr->prev != NULL) 
    {
      n->attr = n->attr->prev;
    }

  s = To_Lower_String(attr_name);
  if(!strcmp(s,"id"))
    {
      XML_Set_Node_Id(n,attr_value);
      printf("\n. Node '%s' id is '%s'",n->name,n->id);
    }
  Free(s);

  return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_attr *XML_Make_Attribute(xml_attr *prev, char *attr_name, char *attr_value)
{
  xml_attr *new_attr;

  new_attr = (xml_attr *)mCalloc(1,sizeof(xml_attr));
  
  new_attr->prev = prev;
  new_attr->next = NULL;
  if(prev != NULL) prev->next = new_attr;
  
  new_attr->name = (char *)mCalloc(strlen(attr_name)+1,sizeof(char));
  strcpy(new_attr->name,attr_name);

  new_attr->value = (char *)mCalloc(strlen(attr_value)+1,sizeof(char));
  strcpy(new_attr->value,attr_value);

  return new_attr;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Set_Node_Id(xml_node *n, char *id)
{
  XML_Make_Node_Id(n,id);
  strcpy(n->id,id);
  return(0);
}
 
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Make_Node_Id(xml_node *n, char *id)
{
  n->id = (char *)mCalloc(strlen(id)+1,sizeof(char));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int XML_Set_Node_Value(xml_node *n, char *val)
{
  XML_Make_Node_Value(n,val);
  strcpy(n->value,val);
  return(0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void XML_Make_Node_Value(xml_node *n, char *val)
{
  n->value = (char *)mCalloc(strlen(val)+1,sizeof(char));
}
