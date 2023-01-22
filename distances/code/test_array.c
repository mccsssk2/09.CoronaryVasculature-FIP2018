#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <string.h>


// Driver method
int main()
{
    int *parent; 
	parent = (int *) calloc(7, sizeof(int));
    int n = sizeof parent / sizeof parent[0];
printf("%d\n", n);
//    Node *root = createTree(parent, n);
//    cout << "Inorder Traversal of constructed tree\n";
//    inorder(root);
//    newLine();
}
