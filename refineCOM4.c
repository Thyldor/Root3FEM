#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"


struct triangle;

struct vertex{
    int tag;
    double x;
    double y;
    int isBorder; //0 if it isn't, 1 if it is
};
struct edge{
    struct vertex *a;
    struct vertex *b;
    struct triangle *triangle1;
    struct triangle *triangle2;
    int isBorder; // 0 false, 1 true
    int isOdd; //0 false; 1 true

};
struct triangle{ //maybe change to arrays of size 3 for the linked entities, but leave it as is for testing
    struct vertex *vertex1;
    struct vertex *vertex2;
    struct vertex *vertex3;
    struct edge *edge1;
    struct edge *edge2;
    struct edge *edge3;
    struct triangle *neighbour1;
    struct triangle *neighbour2;
    struct triangle *neighbour3;
    int type;
};



int edgeIndex = 0, triangleIndex = 0, vertexIndex = 0;
struct vertex *vertices;
struct edge *edges;
struct triangle *triangles;


void refineTriangle(struct triangle * t);
void evenRefine(struct triangle * t);



struct edge * findEdge(int i1, int i2){
    for(int i=0; i<edgeIndex; i++){
        if(((*edges[i].a).tag == i1 && (*edges[i].b).tag == i2)  ||  ((*edges[i].b).tag == i1 && (*edges[i].a).tag == i2))
            return &edges[i];
    }
    return NULL;
}

struct vertex * findThirdVertex(struct triangle * t, struct edge *e){
    if(t->vertex1->tag != e->a->tag && t->vertex1->tag != e->b->tag) return t->vertex1;
    if(t->vertex2->tag != e->a->tag && t->vertex2->tag != e->b->tag) return t->vertex2;
    return t->vertex3;
}
void linkTriangles(struct triangle* t1, struct triangle* t2, struct edge* e){//we need to know the edge we're linking them along
    if(t1->edge1 == e) t1->neighbour1= t2;
    if(t1->edge2 == e) t1->neighbour2= t2;
    if(t1->edge3 == e) t1->neighbour3= t2;
    if(t2->edge1 == e) t2->neighbour1= t1;
    if(t2->edge2 == e) t2->neighbour2= t1;
    if(t2->edge3 == e) t2->neighbour3= t1;
}
void updateEdgeNeighbour(struct edge * e, struct triangle * newT, struct triangle * oldT){
    if(e->triangle1 == oldT){
        e->triangle1 = newT;
        return;
    }
    if(e->triangle2 == oldT){
        e->triangle2 = newT;
    }
}
struct edge * findEdgeInTriangle(struct triangle *t, struct vertex *v1, struct vertex *v2){
    if((t->edge1->a == v1 || t->edge1->b == v1)&&(t->edge1->a == v2 || t->edge1->b == v2))
            return t->edge1;
    if((t->edge2->a == v1 || t->edge2->b == v1)&&(t->edge2->a == v2 || t->edge2->b == v2))
        return t->edge2;
    return t->edge3;
}
void updateTriangleEdge(struct triangle *t, struct edge *newE, struct edge *oldE){
    if(t->edge1 == oldE) {
        t->edge1 = newE;
        return;
    }
    if(t->edge2 == oldE) {
        t->edge2 = newE;
        return;
    }
    if(t->edge3 == oldE){
        t->edge3 = newE;
        return;
    }
}
void updateTriangleVertex(struct triangle *t, struct vertex *newV, struct vertex *oldV){
    if(t->vertex1 == oldV){
        t->vertex1 = newV;
        return;
    }
    if(t->vertex2 == oldV){
        t->vertex2 = newV;
        return;
    }
    if(t->vertex3 == oldV){
        t->vertex3 = newV;
        return;
    }
}
void updateTriangleNeighbour(struct triangle *t, struct triangle *newT, struct triangle *oldT){
    if(t->neighbour1 == oldT){
        t->neighbour1 = newT;
        return;
    }
    if(t->neighbour2 == oldT){
        t->neighbour2 = newT;
        return;
    }
    if(t->neighbour3 == oldT){
        t->neighbour3 = newT;
        return;
    }

}
void refreshTriangleNeighbourFromUpdatedEdge(struct edge * e){
    struct triangle *t1 = e->triangle1, *t2 = e->triangle2;


    if(t1!=NULL) {
        if(t1->edge1 == e) t1->neighbour1 = t2;
        else if(t1->edge2 == e) t1->neighbour2 = t2;
        else t1->neighbour3 = t2;
    }

    if(t2!= NULL){
        if(t2->edge1 == e) t2->neighbour1 = t1;
        else if(t2->edge2 == e) t2->neighbour2 = t1;
        else t2->neighbour3 = t1;
    }
}// could be made redundant using linkTriangles
struct vertex * findCommonVertex(struct edge * e1, struct edge * e2){
    if(e1->a->tag == e2->a->tag)
        return e1->a;
    else if(e1->b->tag == e2->a->tag)
        return e1->b;
    return e2->b;
}

struct triangle * findSecondNeighbourAlongEdge(struct edge * e, struct triangle * t){
    if(e->triangle1 == t)
        return e->triangle2;
    return e->triangle1;
}

struct edge * findRefinementEdge(struct triangle *t){
    if(t->type == 1) return t->edge1;
    if(t->type == 2) return t->edge2;
    if(t->type == 3) return t->edge3;
    printf("Why are you trying to find the refinement edge of an even trangle?!\n");
    return NULL;
}

int isObtuse(struct triangle *t){//returns the edge index opposite the obtuse angle if it exists
    double a1, a2, a3;//inner product at each vertex
    
    a2 = (t->vertex1->x - t->vertex2->x)*(t->vertex3->x - t->vertex2->x) + (t->vertex1->y - t->vertex2->y)*(t->vertex3->y - t->vertex2->y);
    if(a2<0) return 3;//this is first, since we expect most obtuse triangles to be of type 3, so we avoid doing pointless calculations
    
    a3 = (t->vertex1->x - t->vertex3->x)*(t->vertex2->x - t->vertex3->x) + (t->vertex1->y - t->vertex3->y)*(t->vertex2->y - t->vertex3->y);
    if(a3<0) return 1;
    
    a1 = (t->vertex2->x - t->vertex1->x)*(t->vertex3->x - t->vertex1->x) + (t->vertex2->y - t->vertex1->y)*(t->vertex3->y - t->vertex1->y); 
    if(a1<0) return 2;
    
    return 0;      
}

void connectCOMS(struct edge * e){
    struct triangle *t1 = e->triangle1, *t2 = e->triangle2;
    struct vertex *com1 = findThirdVertex(t1, e), *com2 = findThirdVertex(t2, e), *v1 = e->a, *v2 = e->b;
    struct edge *e1, *e2;
    
    
    e->a = com1; e->b = com2; e->isOdd = 0;

    e1 = findEdgeInTriangle(t1, com1, v2); e2 = findEdgeInTriangle(t2, com2, v1);
    
    
    //Checking inputs to avoid errors
    if((t1 == NULL) || (t2 == NULL)){
        mexErrMsgIdAndTxt("bitmarker:refineCOM", "connectCOMS: NULL triangle");
    }
    if((com1 == NULL) ||(com2 == NULL)){
        mexErrMsgIdAndTxt("bitmarker:refineCOM", "connectCOMS: NULL third vertices");
    }
    if((e1 == NULL) ||(e2 == NULL)){
        mexErrMsgIdAndTxt("bitmarker:refineCOM", "connectCOMS: NULL connecting edges");
    }   
    
    
    updateTriangleEdge(t1, e2, e1); updateTriangleEdge(t2, e1, e2);
    updateTriangleVertex(t1, com2, v2); updateTriangleVertex(t2, com1, v1);
    updateEdgeNeighbour(e1, t2, t1); updateEdgeNeighbour(e2, t1, t2);
    refreshTriangleNeighbourFromUpdatedEdge(e1);refreshTriangleNeighbourFromUpdatedEdge(e2);
    t1->type = 0; t2->type = 0;
    
    
    //Refine corners again
    //if(t1->edge1->isBorder + t1->edge2->isBorder + t1->edge3->isBorder == 2){
    //    evenRefine(t1);
    //}
    //if(t2->edge1->isBorder + t2->edge2->isBorder + t2->edge3->isBorder == 2){
    //    evenRefine(t2);
    //}
}



void evenRefine(struct triangle *t){
    struct edge *e1 = t->edge1, *e2 = t->edge2, *e3 = t->edge3;
    struct vertex *v1 = findCommonVertex(e1, e3), *v2 = findCommonVertex(e2,e1), *v3 = findCommonVertex(e2, e3);
    struct triangle *t1 = t->neighbour1, *t2 = t->neighbour2, *t3 = t->neighbour3;


    struct vertex *com = &vertices[++vertexIndex];
    struct edge *newE1 = &edges[edgeIndex++], *newE2 = &edges[edgeIndex++], *newE3 = &edges[edgeIndex++];
    struct triangle *newT1 = t, *newT2 = &triangles[triangleIndex++], *newT3 = &triangles[triangleIndex++];


    //Create the centre of mass
    *com = (struct vertex){vertexIndex, (v1->x + v2->x + v3->x)/3, (v1->y + v2->y + v3->y)/3,0};

    //Create the new edges and triangles
    *newE1 = (struct edge) {v1, com, newT1, newT3, 0, 0};
    *newE2 = (struct edge) {v2, com, newT1, newT2, 0, 0};
    *newE3 = (struct edge) {v3, com, newT2, newT3, 0, 0};

    *newT1 = (struct triangle) {v1, com, v2, newE1, newE2, e1, newT3, newT2, t1, 3};
    *newT2 = (struct triangle) {v2, com, v3, newE2, newE3, e2, newT1, newT3, t2, 3};
    *newT3 = (struct triangle) {v3, com, v1, newE3, newE1, e3, newT1, newT2, t3, 3};


    //Updating the old edges (odd and new neighbour)
    //updateEdgeNeighbour(e1, newT1, t); refreshTriangleNeighbourFromUpdatedEdge(e1);//e1's neighbour is already updated to to overwriting t's memory with newT1 since I'm using arrays. This could change if using linked lists
    updateEdgeNeighbour(e2, newT2, t); refreshTriangleNeighbourFromUpdatedEdge(e2);
    updateEdgeNeighbour(e3, newT3, t); refreshTriangleNeighbourFromUpdatedEdge(e3);


    
    if(isObtuse(newT1) == 3){  //checking if the newly created triangles are obtuse  
        if(e1->isOdd) connectCOMS(e1); //if it is odd, that means that the neighbouring triangle has already been refined, so we would need to connect the two COMs
        else e1->isOdd = 1;
    }else newT1->type = 0;

    if(isObtuse(newT2) == 3){ 
        if(e2->isOdd) connectCOMS(e2);
        else e2->isOdd = 1;
    }else newT2->type = 0;

    if(isObtuse(newT3) == 3){ 
        if(e3->isOdd) connectCOMS(e3);
        else e3->isOdd = 1;
    }else newT3->type = 0;
}


void oddRefine(struct triangle *t){
    struct edge * refinementEdge = findRefinementEdge(t);
    int obt;
    
    
    if (refinementEdge == NULL){
          mexErrMsgIdAndTxt("bitmarker:refineCOM", "Null refinement edge in odd refinement");
    }
    
    if(refinementEdge->isBorder){
        struct edge *e1 ,*e2 , *e3 = refinementEdge;
        if(t->type != 1){
            e1 = t->edge1;
            if(t->type != 2)
                e2 = t->edge2;
            else e2 = t->edge3;
        }else{
            e1 = t->edge2;
            e2 = t->edge3;
        }
        
        struct vertex *v1 = findCommonVertex(e1, e3), *v2 = findCommonVertex(e1, e2), *v3 = findCommonVertex(e2, e3);
        struct triangle *e1Neighbour = findSecondNeighbourAlongEdge(e1, t), *e2Neighbour = findSecondNeighbourAlongEdge(e2, t);

        
        struct vertex * newV1 = &vertices[++vertexIndex], *newV2 = &vertices[++vertexIndex];
        struct edge *newE1 = &edges[edgeIndex++], *newE2 = refinementEdge, *newE3 = &edges[edgeIndex++], *newE4 = &edges[edgeIndex++], *newE5 = &edges[edgeIndex++];
        struct triangle *newT1 = &triangles[triangleIndex++], *newT2 = t, *newT3 = &triangles[triangleIndex++];


        *newV1 = (struct vertex){vertexIndex-1, (2*(v1->x) + v3->x)/3, (2*(v1->y) + v3->y)/3, 1}; //creating the first new vertex at one third of the border edge
        *newV2 = (struct vertex){vertexIndex, (v1->x + 2*(v3->x))/3, (v1->y + 2*(v3->y))/3, 1}; //creating the second new vertex at two thirds of the border edge

        *newE1 = (struct edge){v1, newV1, newT1, NULL, 1, 0};
        *newE2 = (struct edge){newV1, newV2, newT2, NULL, 1, 0};
        *newE3 = (struct edge){newV2, v3, newT3, NULL, 1, 0};
        *newE4 = (struct edge){v2, newV1, newT1, newT2, 0, 0};
        *newE5 = (struct edge){v2, newV2, newT2, newT3, 0, 0};

        *newT1 = (struct triangle){v1, newV1, v2, newE1, newE4, e1, NULL, newT2, e1Neighbour, 3};//the type numbers are the ones expected in the normal case
        *newT2 = (struct triangle){v2, newV1, newV2, newE4, newE2, newE5, newT1, NULL, newT3, 0};//checking for exceptions is done slightly lower
        *newT3 = (struct triangle){v2, newV2, v3, newE5, newE3, e2, newT2, NULL, e2Neighbour, 3};
        

        updateEdgeNeighbour(e1, newT1, t);
        updateEdgeNeighbour(e2, newT3, t);
        updateTriangleNeighbour(e1Neighbour, newT1, t);
        updateTriangleNeighbour(e2Neighbour, newT3, t);
        
        
        

        //check if edges are odd and merge if yes
        if(isObtuse(newT1) == 3){ 
            if(e1->isOdd) connectCOMS(e1);
            else e1->isOdd = 1;
        }else newT1->type = 0;//this is the exceptional case
        
        if(isObtuse(newT3) == 3){ 
        	if(e2->isOdd) connectCOMS(e2);
            else e2->isOdd = 1;
        }else newT3->type = 0;//this is the exceptional case 
        
        obt = isObtuse(newT2);
        if(obt!=0){ //this is the exceptional case
            if(obt == 1)newE4->isOdd = 1;
            else if(obt == 2)newE2->isOdd = 1;
            else newE5->isOdd = 1;
        }
        
    }
    
    
    else{//if the triangle is an interior obtuse triangle
        struct triangle * refinementNeighbour;
                
        while(t->type != 0){//if the triangle is still odd, it means its neighbour was also odd, so one refinement would only make it even. A second refinement is needed to refine the current triangle
            refinementNeighbour = findSecondNeighbourAlongEdge(refinementEdge, t);
            if(refinementNeighbour == NULL){
                mexErrMsgIdAndTxt("bitmarker:refineCOM", "Couldn't find refinement neighbour for the second time");
            }
            refineTriangle(refinementNeighbour);
        }
    }
}


void refineTriangle(struct triangle *t){

    if(t->type == 0){//even type/ type A/ type 0
        evenRefine(t);
    }else{//odd type/ type B / type 1
        oddRefine(t);
    }
}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    mxDouble *mxVertices, *mxTriangles, *mxBorderEdges, *mxMarkedTriangles, *mxTriangleType;
    int i,j,k,i1,i2,i3, bEdges, noOfTs,v1,v2, typeIndex, tp;
    double x1,x2,x3,y1,y2,y3,m;
    
    if(nrhs !=6){
        mexErrMsgIdAndTxt("bitmarker:refineCOM", "There should be 6 inputs");
    }    
    if((mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) || (mxGetN(prhs[0]) != 2)){
        mexErrMsgIdAndTxt("bitmarker:refineCOM", "Wrong vertex matrix input (input 1). It should be an Nx2 double matrix");
    }
    if((mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) || (mxGetN(prhs[1]) != 3)){
        mexErrMsgIdAndTxt("bitmarker:refineCOM", "Wrong triangle matrix input (input 2). It should be an Nx3 double matrix");
    }
    if((mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) || (mxGetN(prhs[2])!=2)){
        mexErrMsgIdAndTxt("bitmarker:refineCOM", "Wrong dirichlet matrix input (input 3). It should be an Nx2 double matrix");
    }
    if((mxGetClassID(prhs[4]) != mxDOUBLE_CLASS) || (mxGetN(prhs[4])!=1)){
        mexErrMsgIdAndTxt("bitmarker:refineCOM", "Wrong marked matrix input (input 5). It should be an Nx1 double matrix");
    }
    if((mxGetClassID(prhs[5]) != mxDOUBLE_CLASS) || (mxGetN(prhs[5])!=1)){
        mexErrMsgIdAndTxt("bitmarker:refineCOM", "Wrong type matrix input (input 6). It should be an Nx1 double matrix");
    }
    
    
    
    vertexIndex = mxGetM(prhs[0]);
    triangleIndex = mxGetM(prhs[1]);
    edgeIndex = 0;//Apparently Matlab keeps global variables in memory, so unless I reset it, it runs out of bounds
    bEdges = mxGetM(prhs[2]);
    noOfTs = mxGetM(prhs[4]);
    typeIndex = mxGetM(prhs[5]);
    
    vertices = malloc(sizeof(*vertices)*(2*vertexIndex + (2*triangleIndex) + 1));//(vertexIndex+2*triangleIndex+1)
    if(vertices == NULL) mexErrMsgIdAndTxt("bitmarker:refineCOM", "Failed vertex memory allocation");
    triangles = malloc(sizeof(*triangles) * 4 * triangleIndex);//3 * triangleIndex
    if(vertices == NULL) mexErrMsgIdAndTxt("bitmarker:refineCOM", "Failed triangle memory allocation");
    edges = malloc(sizeof(*edges) * 4 * ((3*triangleIndex + bEdges)/2)); //3 * ((3*triangleIndex + bEdges)/2)
    if(vertices == NULL) mexErrMsgIdAndTxt("bitmarker:refineCOM", "Failed edge memory allocation");
    
    mxVertices = mxGetDoubles(prhs[0]);
    mxTriangles = mxGetDoubles(prhs[1]);
    mxBorderEdges = mxGetDoubles(prhs[2]);
    mxMarkedTriangles = mxGetDoubles(prhs[4]);
    mxTriangleType = mxGetDoubles(prhs[5]);
    
    for(i = 1; i<=vertexIndex; i++){
        vertices[i] = (struct vertex){i, mxVertices[i-1], mxVertices[i-1+vertexIndex], 0};
    }
    
    for(i=0; i<triangleIndex; i++){
        i1 = (int) mxTriangles[i];
        i2 = (int) mxTriangles[triangleIndex+i];
        i3 = (int) mxTriangles[2*triangleIndex+i];
        triangles[i] = (struct triangle) {&vertices[i1], &vertices[i2], &vertices[i3], findEdge(i1, i2), findEdge(i2, i3),
                        findEdge(i1, i3), NULL, NULL, NULL, 0}; //we'll check for neighbours later, since we need to check the edges first



        if(triangles[i].edge1 == NULL){ //Edge didn't already exist, so we make a new one. This also means its neighbour hasn't been created yet
            edges[edgeIndex] = (struct edge) {&vertices[i1], &vertices[i2], &triangles[i], NULL, 0, 0};
            triangles[i].edge1 = &edges[edgeIndex];
            edgeIndex++;
        }else{//If the edge already existed, that means its neighbour also already exists, so we need to link the two
            triangles[i].edge1->triangle2 = &triangles[i];
            linkTriangles((*triangles[i].edge1).triangle1, &triangles[i], triangles[i].edge1);
        }
        if(triangles[i].edge2 == NULL){
            edges[edgeIndex] = (struct edge) {&vertices[i2], &vertices[i3], &triangles[i], NULL, 0, 0};
            triangles[i].edge2 = &edges[edgeIndex];
            edgeIndex++;
        }else{
            (*triangles[i].edge2).triangle2 = &triangles[i];
            linkTriangles((*triangles[i].edge2).triangle1, &triangles[i], triangles[i].edge2);
        }
        if(triangles[i].edge3 == NULL){
            edges[edgeIndex] = (struct edge){&vertices[i1], &vertices[i3], &triangles[i], NULL, 0, 0};
            triangles[i].edge3 = &edges[edgeIndex];
            edgeIndex++;
        }else{
            triangles[i].edge3->triangle2 = &triangles[i];
            linkTriangles(triangles[i].edge3->triangle1, &triangles[i], triangles[i].edge3);
        }
    } 
    
    
    //triangle type initialization
    if(typeIndex == triangleIndex)//triangle type matrix provided
        for(i=0; i<triangleIndex; i++){
            triangles[i].type = (int)mxTriangleType[i];
            if(triangles[i].type != 0){
                struct edge * refinementEdge = findRefinementEdge(&triangles[i]);
                if(refinementEdge->isOdd == 0){
                    refinementEdge->isOdd = 1;
                }
                else printf("This should never happen given a nice mesh\n");//could add a connectCOMS case, to avoid issues with initial meshes with opposing obtuse triangles
            }
        }
    else{ //triangle type matrix not provided
        for(i=0; i<triangleIndex; i++){//could add a connectCOMS case, to avoid issues with initial meshes with opposing obtuse triangles
            tp = isObtuse(&triangles[i]);
            triangles[i].type = tp;
            if(tp == 1) triangles[i].edge1->isOdd = 1;
            else if(tp == 2) triangles[i].edge2->isOdd = 1;
            else if(tp == 3) triangles[i].edge3->isOdd = 1;            
        }
    }
    
    
    for(i=0;i<bEdges;i++){
        i1 = (int) mxBorderEdges[i];
        i2 = (int) mxBorderEdges[bEdges+i];
        if(findEdge(i1,i2)!=NULL){//redundancy to not cause crashes
            findEdge(i1,i2)->isBorder = 1;
            vertices[i1].isBorder = 1;
            vertices[i2].isBorder = 1;
        }
    }
    
    //if((triangleType == NULL) || (triangleIndex == 8) || (triangleIndex == 25)){//Kind of a hack
    //    mexErrMsgIdAndTxt("bitmarker:refineCOM", "This should never happen if inputting triangly type param\n");
    //    free(triangleType);
    //    triangleType = calloc(triangleIndex, sizeof(int));
    //}
    
    
    
    
    //Go through the marked list and refine all those triangles once, avoiding double refinements
    int marked[3*noOfTs];
    
    for(i=0; i< noOfTs; i++){
        marked[3*i] = triangles[((int)mxMarkedTriangles[i])-1].vertex1->tag;
        marked[3*i+1] = triangles[((int)mxMarkedTriangles[i])-1].vertex2->tag;
        marked[3*i+2] = triangles[((int)mxMarkedTriangles[i])-1].vertex3->tag;
    }
    
    for(i=0; i< noOfTs; i++){
        if((triangles[((int)mxMarkedTriangles[i])-1].vertex1->tag == marked[3*i]) &&
           (triangles[((int)mxMarkedTriangles[i])-1].vertex2->tag == marked[3*i+1]) &&
           (triangles[((int)mxMarkedTriangles[i])-1].vertex3->tag == marked[3*i+2]))
                refineTriangle(&triangles[((int)mxMarkedTriangles[i])-1]);
    }
    
    
    
    
    //---------------------------------------------------------------------
    // ----------------------------OUTPUT----------------------------------
    //---------------------------------------------------------------------
    
    plhs[3] = mxCreateDoubleMatrix(0,0,0);
    
    //-------------------Setting border edge output------------------------
    bEdges = 0;
    for(i=0; i<edgeIndex; i++){//could be replaced with a global border edge count instead of looping over all edges
        if(edges[i].isBorder) bEdges++;
    }
    int bEdgesTags[3][bEdges], end1=0, end2=0;//bEdgesTag = [tag1, tag2, hasBeenUsed(T/F)]
    j=0;
    for(i=0; i<edgeIndex; i++){
        if(edges[i].isBorder){ 
            bEdgesTags[0][j]=edges[i].a->tag;
            bEdgesTags[1][j]=edges[i].b->tag;
            bEdgesTags[2][j++] = 0;            
        }
    }
    
    plhs[2] = mxCreateDoubleMatrix(bEdges, 2, mxREAL);
    mxBorderEdges = mxGetDoubles(plhs[2]);
    
    j=0;i=0;
    while(j<bEdges){
        if(end1 == end2){//ends should be equal only at the start or if there is more than one loop of edges
            while(bEdgesTags[2][i]){//find an unused edge (edge that isn't in one of the alreadt existin loops)
                i=(i+1)%bEdges;
            }
            end1 = bEdgesTags[0][i];
            end2 = bEdgesTags[1][i];
            mxBorderEdges[j] = bEdgesTags[0][i];
            mxBorderEdges[bEdges + j] = bEdgesTags[1][i];
            bEdgesTags[2][i] = 1;
            j++; i=(i+1)%bEdges;
            continue;
        }
        
        if(bEdgesTags[0][i] == end1){
            end1 = bEdgesTags[1][i];
            mxBorderEdges[j] = bEdgesTags[1][i];
            mxBorderEdges[bEdges + j] = bEdgesTags[0][i];
            bEdgesTags[2][i] = 1;
            j++; i=(i+1)%bEdges;
            continue;
        }
        
        if(bEdgesTags[1][i] == end1){
            end1 = bEdgesTags[0][i];
            mxBorderEdges[j] = bEdgesTags[0][i];
            mxBorderEdges[bEdges + j] = bEdgesTags[1][i];
            bEdgesTags[2][i] = 1;
            j++; i=(i+1)%bEdges;
            continue;
        }
        
        if(bEdgesTags[0][i] == end2){
            end2 = bEdgesTags[1][i];
            mxBorderEdges[j] = bEdgesTags[0][i];
            mxBorderEdges[bEdges + j] = bEdgesTags[1][i];
            bEdgesTags[2][i] = 1;
            j++; i=(i+1)%bEdges;
            continue;
        }
        
        if(bEdgesTags[1][i] == end2){
            end2 = bEdgesTags[0][i];
            mxBorderEdges[j] = bEdgesTags[1][i];
            mxBorderEdges[bEdges + j] = bEdgesTags[0][i];
            bEdgesTags[2][i] = 1;
            j++; i=(i+1)%bEdges;
            continue;
        }
        
        i = (i+1)%bEdges;
    }
    
    


    //---------------------Setting triangle output-------------------------
   
    plhs[4] = mxCreateDoubleMatrix(triangleIndex, 1, mxREAL);
    mxTriangleType = mxGetDoubles(plhs[4]);
    
    plhs[1] = mxCreateDoubleMatrix(triangleIndex, 3, mxREAL);
    mxTriangles = mxGetDoubles(plhs[1]);
    struct edge * refinementEdge;
    int tag1,tag2;
    
    for(i=0; i<triangleIndex; i++){
        x1 = triangles[i].vertex1->x; y1 = triangles[i].vertex1->y;
        x2 = triangles[i].vertex2->x; y2 = triangles[i].vertex2->y;
        x3 = triangles[i].vertex3->x; y3 = triangles[i].vertex3->y;
        
        
        
        if(x1 == x2){
            if(x3<x1){
                if(y1>y2){
                    mxTriangles[i] = triangles[i].vertex1->tag;
                    mxTriangles[triangleIndex+i] = triangles[i].vertex3->tag;
                    mxTriangles[2*triangleIndex+i] = triangles[i].vertex2->tag;
                }
                else{
                    mxTriangles[i] = triangles[i].vertex1->tag;
                    mxTriangles[triangleIndex+i] = triangles[i].vertex2->tag;
                    mxTriangles[2*triangleIndex+i] = triangles[i].vertex3->tag;
                }
            }
            else{
                if(y1>y2){
                    mxTriangles[i] = triangles[i].vertex1->tag;
                    mxTriangles[triangleIndex+i] = triangles[i].vertex2->tag;
                    mxTriangles[2*triangleIndex+i] = triangles[i].vertex3->tag;
                }
                else{
                    mxTriangles[i] = triangles[i].vertex1->tag;
                    mxTriangles[triangleIndex+i] = triangles[i].vertex3->tag;
                    mxTriangles[2*triangleIndex+i] = triangles[i].vertex2->tag;                    
                }
            }
        }
        
        
        
        else{
        
        m = (y2-y1)/(x2-x1);
        if(y3 > m*x3 + y1 - (m * x1)){
            if(x1 < x2){
                mxTriangles[i] = triangles[i].vertex1->tag;
                mxTriangles[triangleIndex+i] = triangles[i].vertex2->tag;
                mxTriangles[2*triangleIndex+i] = triangles[i].vertex3->tag;
            }
            else{
                mxTriangles[i] = triangles[i].vertex1->tag;
                mxTriangles[triangleIndex+i] = triangles[i].vertex3->tag;
                mxTriangles[2*triangleIndex+i] = triangles[i].vertex2->tag;
            }
        }
        else{
            if(x1 < x2){
                mxTriangles[i] = triangles[i].vertex1->tag;
                mxTriangles[triangleIndex+i] = triangles[i].vertex3->tag;
                mxTriangles[2*triangleIndex+i] = triangles[i].vertex2->tag;
            }
            else{
                mxTriangles[i] = triangles[i].vertex1->tag;
                mxTriangles[triangleIndex+i] = triangles[i].vertex2->tag;
                mxTriangles[2*triangleIndex+i] = triangles[i].vertex3->tag;
            }
        }
        }
        
        if(triangles[i].type){
            refinementEdge = findRefinementEdge(&triangles[i]);
            tag1 = refinementEdge->a->tag;
            tag2 = refinementEdge->b->tag;
            
            if(mxTriangles[i] == tag1 || mxTriangles[i] == tag2){
                if(mxTriangles[triangleIndex + i] == tag1 || mxTriangles[triangleIndex + i] == tag2){
                    mxTriangleType[i] = 1;
                }
                else mxTriangleType[i] = 3;
            }
            else mxTriangleType[i] = 2;
        }
        else mxTriangleType[i] = 0;
        
    }
    
    free(edges);
    free(triangles);
    
    //----------------------Setting vertex output--------------------------
    plhs[0] = mxCreateDoubleMatrix(vertexIndex, 2, mxREAL);
    mxVertices = mxGetDoubles(plhs[0]);
    
    for(i=0; i<vertexIndex; i++){
        mxVertices[i] = vertices[i+1].x;
        mxVertices[vertexIndex + i] = vertices[i+1].y;
    }  
    free(vertices);
    
}






