#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct node {
    struct node *left;
    struct node *right;
    char *string;
} node;

node *root; /* automatically initialized to NULL */

node *make_node(char const *string) {
    node *ret = malloc(sizeof(node));
    if (ret == NULL)
        return NULL;
    ret->string = malloc(strlen(string) + 1);
    if (ret->string == NULL) {
        free(ret);
        return NULL;
    }
    strcpy(ret->string, string);
    ret->left = NULL;
    ret->right = NULL;
    return ret;
}

void del_node(node *node) {
    free(node->string);
    free(node);
}

int insert(const char *string, node **root) {
    if (*root == NULL)
        *root = make_node(string);
    else {
        node *iter = *root;
        for(;;) {
            int cmp = strcmp(string, iter->string);
            if ( 0 == cmp)
                /* duplicate string - ignore it. */
                return 0;
            else if (1 == cmp) {
                if (NULL == iter->right ) {
                    iter->right = make_node(string);
                    break;
                }
                else iter=iter->right;
            }
            else if ( NULL == iter->left ) {
                iter->left = make_node(string);
                break;
            }
            else
                iter = iter->left;
        }
    }
    return 1;
}

void print(node *root) {
    if (root->left != NULL )
        print(root->left);
    fputs(root->string, stdout);
    if ( root->right != NULL )
        print(root->right);
}

int main() {
    char line[100];

    while (fgets(line, 100, stdin))
        insert(line, &root);
    print(root);
    return 0;
}
