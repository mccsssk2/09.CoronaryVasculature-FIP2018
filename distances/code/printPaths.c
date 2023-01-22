    /*
     * C Program to Print all the Paths from the Root to the Leaf 
     * in a Tree 
Output:
 cc tree20.c
$ a.out
40 -> 20 -> 10 ->
40 -> 20 -> 30 ->
40 -> 60 -> 80 -> 90 ->

40 -> 20 -> 10 -> 
40 -> 20 -> 30 -> 
40 -> 60 -> 80 -> 90 -> 

     */

    #include<stdio.h>
    #include<stdlib.h>
    #include <string.h>

    struct node{
       int data;
       struct node* left;
       struct node* right;
  };

    struct node* newnode(int data){
      struct node* node = (struct node*) malloc(sizeof(struct node));
      node->data = data;
      node->left = NULL;
      node->right = NULL;
      return(node);
    }

    /*_________________________________________________________________________________
____________________________________________________________________________________*/
    
    /*Function which helps the print_path to recursively print all the nodes*/ 
    void print_paths_recur(struct node* node, int path[], int path_len, int n1, int n2){ // this needs to return a value to print paths.
      int i;
      int count;
      int return_value;
      return_value = -1;
      FILE *parentTest;
      if (node == NULL) return; 
      path[path_len] = node->data;
      path_len++;
      if (node->left == NULL && node->right == NULL) {
// this is where you test if n1 and n2 are in the same path or not.	      
//      for (i = 0; i < path_len; i++)   printf("%d -> ", path[i]); printf("\n");    
	      count = 0;
	      for (i = 0; i < path_len; i++){ if(path[i]==n1) count++; if(path[i]==n2) count++; }
              if(count > return_value) return_value = count;
              parentTest = fopen("p.txt","a+"); // temporary file.
              fprintf(parentTest, "%d\n", return_value);
              fclose(parentTest);
      }
      else{
        print_paths_recur(node->left,    path, path_len, n1, n2);    //recursively calls the left node of the tree
        print_paths_recur(node->right, path, path_len, n1, n2);    //recursively calls the right node of the tree
      }
    } // end of print_paths_recur.

    /*Function to store all the paths from the root node to all leaf nodes in  a array*/
    void print_paths(struct node* node, int n1, int n2){  // this needs to return a value to the main. 0 if n1 and n2 are parents
      int path[1000];
      print_paths_recur(node, path, 0, n1, n2);
    }

    /*_________________________________________________________________________________
____________________________________________________________________________________*/
    

    int main(){ 
       /*
        The input tree is as shown below
                    40
                    / \
                20        60
                / \       \
            10        30      80
                              \
                                90
      */ 

int count;

      struct node *root = newnode(40);
      FILE *pTest;
      int count_value, count_max;
      char command[50];
      
      root->left 					= newnode(20);
      root->right 				= newnode(60);
      root->left->left 				= newnode(10);
      root->left->right 			= newnode(30);
      root->right->right 			= newnode(80);
      root->right->right->right 	= newnode(90);


   strcpy( command, "rm -f p.txt" ); // temporary file.
   system(command);
      
   print_paths(root, 20, 30); // this needs a return value so I know in main that the nodes are parents or not. 2 more arguments: int n1, int n2.
 
pTest = fopen("p.txt","r");
count_max = -1;
while(fscanf(pTest,"%d", &count_value)!=EOF) if(count_max < count_value) count_max = count_value;

printf("Parents are 2, non-parents are less than 2: %d\n", count_max);
return 0;
    }
