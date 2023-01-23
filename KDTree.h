#pragma once

using namespace std;
namespace Tmpl8{
class KDTree
{
public: 
    const int k = 3;
    float smallestDist = 1e30f;
    int count = 0;
    // A structure to represent node of kd tree
    struct Node
    {
        int idx = -1;
        float3 point; // To store k dimensional point
        float3 normal;
        Node* left, * right;
    };

    Node* rootNode = nullptr;
    Node* nearestNode = nullptr;
    // A method to create a node of K D tree
    struct Node* newNode(float3 arr, float3 normal)
    {
        struct Node* temp = new Node;

        for (int i = 0; i < k; i++)
            temp->point[i] = arr[i];

        temp->left = temp->right = NULL;
        temp->idx = count;
        temp->normal = normal;
        count++;
        return temp;
    }

    // Inserts a new node and returns root of modified tree
    // The parameter depth is used to decide axis of comparison
    Node* insertRec(Node* root, float3 point, float3 normal, unsigned depth)
    {
        // Tree is empty?
        if (root == NULL){
            Node* temp = newNode(point, normal);
            if (rootNode == nullptr) rootNode = temp;
            return temp;
        }

        // Calculate current dimension (cd) of comparison
        unsigned cd = depth % 3;

        // Compare the new point with root on current dimension 'cd'
        // and decide the left or right subtree
        if (point[cd] < (root->point[cd]))
            root->left = insertRec(root->left, point, normal, depth + 1);
        else
            root->right = insertRec(root->right, point, normal, depth + 1);

        return root;
    }

    // Function to insert a new point with given point in
    // KD Tree and return new root. It mainly uses above recursive
    // function "insertRec()"
    Node* insert(Node* root, float3 point, float3 normal)
    {
        insertRec(root, point, normal, 0);
        return rootNode;
    }

    // A utility method to determine if two Points are same
    // in K Dimensional space
    bool arePointsSame(float3 point1, float3 point2)
    {
        // Compare individual pointinate values
        for (int i = 0; i < k; ++i)
            if (point1[i] != point2[i])
                return false;

        return true;
    }

    float getNearestDist(Node* root, float3 currPoint, int depth) {
        smallestDist = 1e30f;
        if (root == nullptr)
        {
            smallestDist = 1e30f;
            nearestNode = root;
            return smallestDist;
        }
        findNearest(root, currPoint, depth);
        //if (nearestNode == nullptr) { nearestNode = root; smallestDist = 0; }
        return smallestDist;
    }

    void findNearest(Node* root,float3 currPoint, int depth) {
        if (root == nullptr)
            return;
        float dis = length(root->point - currPoint);
        if (nearestNode == nullptr || dis < smallestDist) {
            smallestDist = dis;
            nearestNode = root;
        }

        if (smallestDist == 0)
            return;
        float dx = root->point[depth] - currPoint[depth];
        depth = (depth + 1) % 3;
        findNearest(dx > 0 ? root->left : root->right, currPoint, depth);
        if (dx * dx >= smallestDist)
            return;
        findNearest(dx > 0 ? root->right : root->left, currPoint, depth);
    }
    // Searches a Point represented by "point[]" in the K D tree.
    // The parameter depth is used to determine current axis.
    bool searchRec(Node* root, float3 point, unsigned depth)
    {
        // Base cases
        if (root == NULL)
            return false;
        if (arePointsSame(root->point, point))
            return true;

        // Current dimension is computed using current depth and total
        // dimensions (k)
        unsigned cd = depth % k;

        // Compare point with root with respect to cd (Current dimension)
        if (point[cd] < root->point[cd])
            return searchRec(root->left, point, depth + 1);

        return searchRec(root->right, point, depth + 1);
    }

    // Searches a Point in the K D tree. It mainly uses
    // searchRec()
    bool search(Node* root, float3 point)
    {
        // Pass current depth as 0
        return searchRec(root, point, 0);
    }

};
}
