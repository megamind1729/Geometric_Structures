/* Fortune's algorithm to find the Voronoi diagram of a point set */

#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <cassert>
#include <optional>
#include <fstream>

struct point;
struct node;
struct ptr_to_circle_event;
struct ptr_to_arc;
struct half_edge;

std::multimap<node, ptr_to_circle_event>::iterator find_custom(node);
double distance(point, point);
std::optional<double> parabola_intersection(point, point, double);
std::pair<half_edge, half_edge> start_edge(point, double, double);
std::optional<point> findCircleEvent(point, point, point);
std::optional<point> findCircumcentre(point, point, point);

void printBeachline();
void printEventQueue();
void printHalfEdges();
void HandleSiteEvent(point);
void HandleCircleEvent(node*, node*);

std::vector<point> input_points;
std::multimap <point, ptr_to_arc> event_queue;
std::multimap <node, ptr_to_circle_event> beachline;
std::vector <half_edge> half_edges;

double sweepLine_y;     // global variable : y-coordinate of the sweep line
double sweepLine_y_1;   // global variable : first (maximum) y-coordinate among the sites
double epsilon = 1e-7;  // global variable : small value to handle floating point errors
double epsilon2 = 1e-5; // global variable : small value to handle floating point errors (for circle event)
int circle_event_index = -1;

struct point
{
    int index;  // index of the point, -1 for circle events
    double x, y;
    int point_type; // 0: site event, 1: circle event

    bool operator < (const point &p) const
    {
        return y > p.y || (y == p.y && x < p.x) || (y == p.y && x == p.x && point_type < p.point_type) || (y == p.y && x == p.x && point_type == p.point_type && index > p.index);
    }

    bool operator == (const point &p) const
    {
        return (abs(x - p.x) < epsilon && abs(y - p.y) < epsilon && index == p.index && point_type == p.point_type);
    }
};

struct half_edge{
    std::optional<double> start_point_x, start_point_y;   // Start point of half edge
    std::optional<double> end_point_x, end_point_y;       // End point of half edge
    double u, v;                        // Vector parallel to the half edge

    half_edge(){
        this->start_point_x = std::nullopt; this->start_point_y = std::nullopt;
        this->end_point_x = std::nullopt; this->end_point_y = std::nullopt;
        this->u = 0; this->v = 0;
    }

    void operator = (const half_edge &h){
        this->start_point_x = h.start_point_x; this->start_point_y = h.start_point_y;
        this->end_point_x = h.end_point_x; this->end_point_y = h.end_point_y;
        this->u = h.u; this->v = h.v;
    }

    void print_half_edge(){
        std::cout<< " start : ( "; 
        if(this->start_point_x.has_value()){ std::cout<< this->start_point_x.value() << ","; } else{ std::cout<< "INF ,"; }
        if(this->start_point_y.has_value()){ std::cout<< this->start_point_y.value() << ")";} else{ std::cout<< "INF )";}
        std::cout<< ", end : (";
        if(this->end_point_x.has_value()){ std::cout<< this->end_point_x.value() << ","; } else{ std::cout<< "INF ,"; }
        if(this->end_point_y.has_value()){ std::cout<< this->end_point_y.value() << ")";} else{ std::cout<< "INF )";}
        std::cout<< ", u = " << this->u << ", v = " << this->v << std::endl;
    }
};

struct ptr_to_circle_event{
    point* p_circle_ptr; point* q_circle_ptr;
    half_edge half_edge_traced;

    ptr_to_circle_event(){
        this->p_circle_ptr = nullptr; this->q_circle_ptr = nullptr;
        this->half_edge_traced = half_edge();
    }

    ptr_to_circle_event(point* p_circle_ptr, point* q_circle_ptr, half_edge half_edge_traced){
        this->p_circle_ptr = p_circle_ptr; this->q_circle_ptr = q_circle_ptr;
        this->half_edge_traced = half_edge_traced;
    }
};

struct node
{
    std::optional<point> p_opt, q_opt; 
    bool isEventPoint; // 0: breakpoint, 1: event point

    node(){
        this-> p_opt = std::nullopt; this-> q_opt = std::nullopt;
        this->isEventPoint = 0;
    }

    node(std::optional<point> p, std::optional<point> q, bool isEventPoint){
        this->p_opt = p; this->q_opt = q;
        this->isEventPoint = isEventPoint;
    }

    bool operator ==(const node &n) const
    {
        // assert(p_opt.has_value());
        if(!q_opt.has_value()){ return false;}
        return (p_opt.value().index == n.p_opt.value().index && q_opt.value().index == n.q_opt.value().index);
        // if(p_opt.has_value() && n.p_opt.has_value()){
        //     if((abs(p_opt.value().x - n.p_opt.value().x) > epsilon) || (abs(p_opt.value().y - n.p_opt.value().y) > epsilon)){ return false; }
        // }
        // else if(p_opt.has_value() || n.p_opt.has_value()){ return false; }
        // if(q_opt.has_value() && n.q_opt.has_value()){
        //     if((abs(q_opt.value().x - n.q_opt.value().x) > epsilon) || (abs(q_opt.value().y - n.q_opt.value().y) > epsilon)){ return false; }
        // }
        // else if(q_opt.has_value() || n.q_opt.has_value()){ return false; }
        // return isEventPoint == n.isEventPoint;
    }

    bool operator < (const node &n) const
    {
        // For debugging purposes,
        // std::cout<< "\nComparing nodes :\n";
        // std::cout<< "This: " << p_opt.has_value() << " " << q_opt.has_value() << " " << isEventPoint << std::endl;
        // if(p_opt.has_value()) std::cout<< "p_opt: (" << p_opt.value().x << ", " << p_opt.value().y << ") ";
        // else std::cout<< "p_opt: NULL ";
        // if(q_opt.has_value()) std::cout<< "q_opt: (" << q_opt.value().x << ", " << q_opt.value().y << ") " << std::endl;
        // else std::cout<< "q_opt: NULL " << std::endl;
        // std::cout<< "That: " << n.p_opt.has_value() << " " << n.q_opt.has_value() << " " << n.isEventPoint << std::endl;
        // if(n.p_opt.has_value()) std::cout<< "p_opt: (" << n.p_opt.value().x << ", " << n.p_opt.value().y << ") ";
        // else std::cout<< "p_opt: NULL ";
        // if(n.q_opt.has_value()) std::cout<< "q_opt: (" << n.q_opt.value().x << ", " << n.q_opt.value().y << ") " << std::endl;
        // else std::cout<< "q_opt: NULL " << std::endl;

        // // assertions to ensure that both p_opt and q_opt are not empty at the same time for any node
        // assert(p_opt.has_value() || q_opt.has_value());
        // assert(n.p_opt.has_value() || n.q_opt.has_value());
        bool ans;

        // If both are breakpoints (Type 1),
        if(!isEventPoint && !n.isEventPoint){
            if(!p_opt.has_value()) {
                if(!n.p_opt.has_value()){ // Both are of the form (NULL, q)
                    // assert(q_opt.has_value() && n.q_opt.has_value());
                    ans = (q_opt.value().x < n.q_opt.value().x);
                    // std::cout<< "Ans1: " << ans << std::endl;
                    return ans;
                }
                else{
                    ans = true;
                    // std::cout<< "Ans2: " << ans << std::endl;
                    return true;
                }
            }
            else if(!q_opt.has_value()){
                if(!n.q_opt.has_value()){ // Both are of the form (p, NULL)
                    // assert(p_opt.has_value() && n.p_opt.has_value());
                    ans = (p_opt.value().x < n.p_opt.value().x);
                    // std::cout<< "Ans3: " << ans << std::endl;
                    return ans;
                }
                else{
                    ans = false;
                    // std::cout<< "Ans4: " << ans << std::endl;
                    return false;
                }
            }
            else if(!n.p_opt.has_value()){ // n is of the form (NULL, q)
                ans = false;
                // std::cout<< "Ans5: " << ans << std::endl;
                return false;
            }
            else if(!n.q_opt.has_value()){ // n is of the form (p, NULL)
                ans = true;
                // std::cout<< "Ans6: " << ans << std::endl;
                return true;
            }            
            else{
                // assert(p_opt.has_value() && q_opt.has_value() && n.p_opt.has_value() && n.q_opt.has_value());
                // std::cout<< "p_opt: (" << p_opt.value().x << ", " << p_opt.value().y << ") " << "q_opt: (" << q_opt.value().x << ", " << q_opt.value().y << ") " << "n.p_opt: (" << n.p_opt.value().x << ", " << n.p_opt.value().y << ") " << "n.q_opt: (" << n.q_opt.value().x << ", " << n.q_opt.value().y << ") " << std::endl;
                std::optional<double> x1_opt = parabola_intersection(p_opt.value(), q_opt.value(), sweepLine_y+epsilon2);
                std::optional<double> x2_opt = parabola_intersection(n.p_opt.value(), n.q_opt.value(), sweepLine_y+epsilon2);
                // assert(x1_opt.value() && x2_opt.value());
                ans = x1_opt.value() < x2_opt.value();
                // std::cout<< "Case 1: Both breakpoints\n";
                // std::cout<< "Ans7: " << ans << std::endl;
                return ans;
            }
        }
        // If this is an event point (Type 2 or 3) and n is a breakpoint (Type 1),
        else if(isEventPoint && !n.isEventPoint){
            // assert(!q_opt.has_value());
            if(!n.p_opt.has_value()){ // n is of the form (NULL, q)
                    ans = false;
                    // std::cout<< "Ans8: " << ans << std::endl;
                    return false;
                }
            else if(!n.q_opt.has_value()){ // n is of the form (p, NULL)
                ans = true;
                // std::cout<< "Ans9: " << ans << std::endl;
                return true;
            }
            else{    // Type 2 and Type 1
                // assert(p_opt.has_value() && n.p_opt.has_value() && n.q_opt.has_value());
                ans = (p_opt.value().x < parabola_intersection(n.p_opt.value(), n.q_opt.value(), sweepLine_y));
                // std::cout << ("Case 2: Event and breakpoint\n");
                // std::cout<< "Ans10: " << ans << std::endl;
                return ans;
            }
        }
        // If this is a breakpoint(Type 1) and n is an event point (Type 2 or 3),
        else if(!isEventPoint && n.isEventPoint){
            // assert(!n.q_opt.has_value());
            if(!p_opt.has_value()){ // This is of the form (NULL, q)
                // std::cout<< "Ans19: " << 1 << std::endl;
                return true;
            }
            else if(!q_opt.has_value()){ // This is of the form (p, NULL)
                // std::cout<< "Ans20: " << 0 << std::endl;
                return false;
            }
            else{
                // assert(p_opt.has_value() && q_opt.has_value() && n.p_opt.has_value());
                ans = (parabola_intersection(p_opt.value(), q_opt.value(), sweepLine_y) < n.p_opt.value().x);
                // std::cout << "Case 3: Breakpoint and event\n";
                // std::cout<< "Ans21: " << ans << std::endl;
                return ans;
            }
        }
        // If both are event points,
        else{
            std::cout<< "HOW? WHY?" << std::endl;
            // assert(!q_opt.has_value() && !n.q_opt.has_value());
            // assert(p_opt.has_value() && n.p_opt.has_value());
            return p_opt.value().x < n.p_opt.value().x;
        }
    }
};

struct ptr_to_arc{
    node* left_end;     // Pointer to the left end of the arc
    node* right_end;    // Pointer to the right end of the arc

    ptr_to_arc(){
        this->left_end = nullptr;
        this->right_end = nullptr;
    }
    ptr_to_arc(node* left_end, node* right_end){
        this->left_end = left_end;
        this->right_end = right_end;
    }
};

void writeToFile(const std::vector<half_edge>& half_edges) {
    std::ofstream file("edge.txt");

    // Write half edges to file
    for (const auto& half_edge : half_edges) {
        file << (half_edge.start_point_x.has_value() ? std::to_string(half_edge.start_point_x.value()) : "null") << " "
             << (half_edge.start_point_y.has_value() ? std::to_string(half_edge.start_point_y.value()) : "null") << " "
             << (half_edge.end_point_x.has_value() ? std::to_string(half_edge.end_point_x.value()) : "null") << " "
             << (half_edge.end_point_y.has_value() ? std::to_string(half_edge.end_point_y.value()) : "null") << " "
             << half_edge.u << " " << half_edge.v << "\n";
    }

    file.close();
}

// Function to find the iterator corresponding to a node in the beachline (called only on nodes whihc exist)
std::multimap<node, ptr_to_circle_event>::iterator find_custom(node n){
    auto it = beachline.upper_bound(n);
    // if(it == beachline.end()){ std::cout<< " Not intended\n"; it--; }
    // std::cout<< it->first.p_opt.has_value() << " " << it->first.q_opt.has_value() << " " << it->first.isEventPoint << "\t";
    // if(it->first.p_opt.has_value()){std::cout<< "p: " << it->first.p_opt.value().index << " (" << it->first.p_opt.value().x << ", " << it->first.p_opt.value().y << ") ";} else{std::cout<< "p: NULL ";}
    // if(it->first.q_opt.has_value()){std::cout<< "q: " << it->first.q_opt.value().index << " (" << it->first.q_opt.value().x << ", " << it->first.q_opt.value().y << ") ";} else{std::cout<< "q: NULL ";}
    // std::cout<< std::endl;
    int i = 2;
    while(it->first.q_opt.has_value() && i > 0){
        it++; i--;
        // std::cout<< it->first.p_opt.has_value() << " " << it->first.q_opt.has_value() << " " << it->first.isEventPoint << "\t";
        // if(it->first.p_opt.has_value()){std::cout<< "p: " << it->first.p_opt.value().index << " (" << it->first.p_opt.value().x << ", " << it->first.p_opt.value().y << ") ";} else{std::cout<< "p: NULL ";}
        // if(it->first.q_opt.has_value()){std::cout<< "q: " << it->first.q_opt.value().index << " (" << it->first.q_opt.value().x << ", " << it->first.q_opt.value().y << ") ";} else{std::cout<< "q: NULL ";}
        // std::cout<< std::endl;
    }
    while(!(it->first == n) && it != beachline.begin()){ 
        it--;
        // std::cout<< it->first.p_opt.has_value() << " " << it->first.q_opt.has_value() << " " << it->first.isEventPoint << "\t";
        // if(it->first.p_opt.has_value()){std::cout<< "p: " << it->first.p_opt.value().index << " (" << it->first.p_opt.value().x << ", " << it->first.p_opt.value().y << ") ";} else{std::cout<< "p: NULL ";}
        // if(it->first.q_opt.has_value()){std::cout<< "q: " << it->first.q_opt.value().index << " (" << it->first.q_opt.value().x << ", " << it->first.q_opt.value().y << ") ";} else{std::cout<< "q: NULL ";} 
        // std::cout<< std::endl;
    }
    // std::cout<< "Find custom over\n";
    return it;
}

// Function to find the distance between two points
// This function will be called to find the circle event
double distance(point p, point q)
{
    return sqrt((p.x - q.x)*(p.x - q.x) + (p.y - q.y)*(p.y - q.y));
}

// Function to find the x-coordinate of the intersection of two parabolas given their foci and the y-coordinate of the directrix
// The intersection point is such that the parabola with focus p is to the left of the parabola with focus q at the intersection point
std::optional<double> parabola_intersection(point p, point q, double l)
{
    double x1 = p.x, y1 = p.y - l, x2 = q.x, y2 = q.y - l;
    // assert(!(x1 == x2 && y1 == y2));

    double a = y2 - y1;
    double b = -2*(x1*y2 - x2*y1);
    double c = y2*x1*x1 - y1*x2*x2 - y1*y2*(y2 - y1);     
    double d = b*b - 4*a*c;

    // For debugging purposes,
    // std::cout<< "\nParabola intersection computation: (Sweep line : "<< l << ")" << std::endl; 
    // std::cout<< "("<< x1 << "," << p.y << ") (" << x2 << "," << q.y << "), " << a << ", " << b << ", " << c << ", " << d << std::endl;

    if(a == 0){ // y1 = y2 and x1 != x2
        if(x1 < x2){
            // std::cout<< "Two parallel parabolas : "<< (x1+x2)/2 << std::endl;
            return (x1+x2)/2;
        }
        else{
            // std::cout<< "No intersection point" << std::endl;
            return std::nullopt;
        }
    }
    
    if(d < epsilon){ d = 0; }
    // assert(d >= 0); // d = 0 is possible when y1 or y2 are zero (that is, sites which are just inserted)
    double ansx_1 = (-b - sqrt(d))/(2*a);
    double ansx_2 = (-b + sqrt(d))/(2*a);
    // std::cout<< "Intersection at : x = " << ansx_2 << std::endl; 
    return ansx_2;
}

// To find the start of the half edges during HandleSiteEvent
std::pair<half_edge, half_edge> start_edge(point p, double x, double l){
    double x1 = p.x, y1 = p.y;
    double y = ((x1-x)*(x1-x) + y1*y1 - l*l)/(2*(y1 - l)); // m = (x-x1)/(y1-l)
    double u_right = y1 - l;
    double v_right = x - x1;
    half_edge h_left = half_edge(); h_left.start_point_x = x; h_left.start_point_y = y; h_left.u = -u_right; h_left.v = -v_right;
    half_edge h_right = half_edge(); h_right.start_point_x = x; h_right.start_point_y = y; h_right.u = u_right; h_right.v = v_right;
    return std::make_pair(h_left, h_right);
}

// To find the circle event (lowermost point of circumcircle of the given points p, q, r)
// Note that the order of arguments cannot be changed
// This function returns the circle event of the arc with left breakpoint (p,q) and right breakpoint (q,r)
std::optional<point> findCircleEvent(point p, point q, point r){
    double x1 = p.x, y1 = p.y, x2 = q.x, y2 = q.y, x3 = r.x, y3 = r.y;

    double B = (x1*x1 + y1*y1)*(y3 - y2) + (x2*x2 + y2*y2)*(y1 - y3) + (x3*x3 + y3*y3)*(y2 - y1);
    double C = (x1*x1 + y1*y1)*(x2 - x3) + (x2*x2 + y2*y2)*(x3 - x1) + (x3*x3 + y3*y3)*(x1 - x2);
    double D = 2*(x1*y2 - x2*y1 + x2*y3 - x3*y2 + x3*y1 - x1*y3);

    // For debugging purposes, 
    // std::cout<< "\nCircle event computation: \n";
    // std::cout<< " ("<< x1 << " " << y1 << ") (" << x2 << "," << y2 << ") (" << x3 << "," << y3 << ") " << B << " " << C << " " << D << std::endl;

    if(D + epsilon >= 0){ // Collinear points or anti-clockwise orientation - No circle event for this arc
        return std::nullopt;
    }
    double x = -B/D;
    double y = -C/D;
    point circumcentre = point({0, x, y, 1});
    double radius = distance(circumcentre, p);
    if(circumcentre.y - radius > sweepLine_y){
        return std::nullopt;
    }
    // std::cout<< "Circle event: (" << x << "," << y << "), Radius = " << radius << std::endl;
    point return_point = point({circle_event_index, x, y-radius, 1});
    circle_event_index--;
    return return_point;
}

// Function to find the circumcentre of focii p, q, r
std::optional<point> findCircumcentre(point p, point q, point r){
    double x1 = p.x, y1 = p.y, x2 = q.x, y2 = q.y, x3 = r.x, y3 = r.y;
    double B = (x1*x1 + y1*y1)*(y3 - y2) + (x2*x2 + y2*y2)*(y1 - y3) + (x3*x3 + y3*y3)*(y2 - y1);
    double C = (x1*x1 + y1*y1)*(x2 - x3) + (x2*x2 + y2*y2)*(x3 - x1) + (x3*x3 + y3*y3)*(x1 - x2);
    double D = 2*(x1*y2 - x2*y1 + x2*y3 - x3*y2 + x3*y1 - x1*y3);
    if(D + epsilon >= 0){
        return std::nullopt;
    }
    double x = -B/D;
    double y = -C/D;
    point circumcentre = point({0, x, y, 1});
    return circumcentre;
}

// Function to display the beachline for debugging purposes
void printBeachline(){
    std::cout << "\n\tBEACHLINE: \n";
    for(auto it = beachline.begin(); it != beachline.end(); it++){
        // assert(it->first.isEventPoint == 0);
        
        std::cout << "(" ;
        if(it->first.p_opt.has_value()){ std::cout << "p: " << it->first.p_opt.value().index << " (" << it->first.p_opt.value().x << ", " << it->first.p_opt.value().y << "), "; }
        else{ std::cout << "p: NULL, "; }
        if(it->first.q_opt.has_value()){ std::cout << "q: " << it->first.q_opt.value().index << " (" << it->first.q_opt.value().x << ", " << it->first.q_opt.value().y << ")) "; }
        else{ std::cout << "q: NULL) "; }
        if(it->second.p_circle_ptr){ std::cout << "(p_circle_ptr: (" << it->second.p_circle_ptr->x << ", " << it->second.p_circle_ptr->y << ") "; }
        else{ std::cout << "(p_circle_ptr: NULL , "; }
        if(it->second.q_circle_ptr){ std::cout << "q_circle_ptr: (" << it->second.q_circle_ptr->x << ", " << it->second.q_circle_ptr->y << ")) "; }
        else{ std::cout << "q_circle_ptr: NULL) "; }
        it->second.half_edge_traced.print_half_edge();
        if(it->first.p_opt.has_value() && it->first.q_opt.has_value()){
            if(parabola_intersection(it->first.p_opt.value(), it->first.q_opt.value(), sweepLine_y).has_value()){
                std::cout << parabola_intersection(it->first.p_opt.value(), it->first.q_opt.value(), sweepLine_y).value() << '\n'<< std::endl;
            }
        }
    }
}

void printEventQueue(){
    std::cout << "\n\tEVENT QUEUE: \n";
    for(auto it = event_queue.begin(); it != event_queue.end(); it++){
        std::cout << it->first.index << " Point: (" << it->first.x << ", " << it->first.y << ") " << " " << it->first.point_type << std::endl;
        if(it->first.point_type == 1){
            // assert(it->second.left_end && it->second.right_end); 
            if(it->second.left_end->p_opt.has_value()){ std::cout << "Left end p : (" << it->second.left_end->p_opt.value().x << "," << it->second.left_end->p_opt.value().y << ") "; }
            else{ std::cout << "Left end p : NULL "; }
            std::cout<< std::endl;
            if(it->second.left_end->q_opt.has_value()){ std::cout << "Left end q : (" << it->second.left_end->q_opt.value().x << "," << it->second.left_end->q_opt.value().y << ") "; }
            else{ std::cout << "Left end q : NULL "; }
            std::cout << std::endl;
            if(it->second.right_end->p_opt.has_value()){ std::cout << "Right end p : (" << it->second.right_end->p_opt.value().x << "," << it->second.right_end->p_opt.value().y << ") "; }
            else{ std::cout << "Right end p : NULL "; }
            std::cout << std::endl;
            if(it->second.right_end->q_opt.has_value()){ std::cout << "Right end q : ("<< it->second.right_end->q_opt.value().x << "," << it->second.right_end->q_opt.value().y << ") "; }
            else{ std::cout << "Right end q : NULL "; }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void printHalfEdges(){
    std::cout << "\n\tHALF EDGES VECTOR: \n";
    for(auto it = half_edges.begin(); it != half_edges.end(); it++){
        std::cout<< " start : ( "; 
        if(it->start_point_x.has_value()){ std::cout<< it->start_point_x.value() << ","; } else{ std::cout<< "INF ,"; }
        if(it->start_point_y.has_value()){ std::cout<< it->start_point_y.value() << ")";} else{ std::cout<< "INF )";}
        std::cout<< ", end : (";
        if(it->end_point_x.has_value()){ std::cout<< it->end_point_x.value() << ","; } else{ std::cout<< "INF ,"; }
        if(it->end_point_y.has_value()){ std::cout<< it->end_point_y.value() << ")";} else{ std::cout<< "INF )";}
        std::cout<< ", u = " << it->u << ", v = " << it->v << std::endl;
    }
}

// Function to handle the site event
void HandleSiteEvent(point p){
    
    // If the beachline is empty, insert the new parabola/focus point into it and return
    if(beachline.empty()){
        std::optional <point> p_opt = p;
        sweepLine_y_1 = sweepLine_y;
        beachline.insert(std::make_pair(node({std::nullopt, p_opt, 0}), ptr_to_circle_event({nullptr, nullptr, half_edge()})));
        beachline.insert(std::make_pair(node({p_opt, std::nullopt, 0}), ptr_to_circle_event({nullptr, nullptr, half_edge()})));
        return;
    }

    // The case where all the parabolas are still degenerate (maximum y among the sites is still same as the current y of the sweep line)
    if(sweepLine_y_1 == sweepLine_y){
        auto it_begin = beachline.begin();
        auto it_end = --beachline.end();
        std::optional<point> q1_opt = it_begin->first.q_opt;
        std::optional<point> q2_opt = it_end->first.p_opt;
        std::optional<point> p_opt = p;
        // assert(q1_opt.has_value() && q2_opt.has_value());
        if(p.x < q1_opt.value().x){
            beachline.erase(it_begin);
            half_edge h = half_edge();
            h.start_point_x = (p.x + q1_opt.value().x)/2; h.start_point_y = std::nullopt; h.u = 0; h.v = -1;
            h.print_half_edge();
            beachline.insert(std::make_pair(node({p_opt, q1_opt, 0}), ptr_to_circle_event({nullptr, nullptr, h})));
            beachline.insert(std::make_pair(node({std::nullopt, p_opt, 0}), ptr_to_circle_event({nullptr, nullptr, half_edge()})));
        }
        else if(p.x > q2_opt.value().x){
            beachline.erase(it_end);
            half_edge h = half_edge();
            h.start_point_x = (p.x + q2_opt.value().x)/2; h.start_point_y = std::nullopt; h.u = 0; h.v = -1;
            std::cout<< "New half edge : \n";
            h.print_half_edge();
            beachline.insert(std::make_pair(node({q2_opt, p_opt, 0}), ptr_to_circle_event({nullptr, nullptr, h})));
            beachline.insert(std::make_pair(node({p_opt, std::nullopt, 0}), ptr_to_circle_event({nullptr, nullptr, half_edge()})));
        }
        else{
            auto it2 = beachline.lower_bound(node({p, std::nullopt, 1}));
            if (it2 != beachline.begin() && it2 != beachline.end()){
                auto it1 = it2; it1--;
                // assert(it1->first.q_opt.has_value() && it2->first.p_opt.has_value());                
                point q2 = it1->first.q_opt.value();
                point q3 = it2->first.p_opt.value();
                // assert(abs(q2.x - q3.x) < epsilon && abs(q2.y - q3.y) < epsilon);
                std::optional<point> q1_opt = it1->first.p_opt;
                std::optional<point> q4_opt = it2->first.q_opt;
                std::optional<point> p_opt = p;
                half_edge h1 = half_edge(); 
                half_edge h2 = half_edge();
                if(p.x > q3.x){
                    std::optional<point> q3_opt = q3;
                    h1.start_point_x = (p.x + q3.x)/2; h1.start_point_y = std::nullopt; h1.u = 0; h1.v = -1;
                    if(q4_opt.has_value()){
                        h2.start_point_x = (p.x + q4_opt.value().x)/2; h2.start_point_y = std::nullopt; h2.u = 0; h2.v = -1;
                    }
                    beachline.erase(it2);
                    beachline.insert(std::make_pair(node({q3_opt, p_opt, 0}), ptr_to_circle_event({nullptr, nullptr, h1})));
                    beachline.insert(std::make_pair(node({p_opt, q4_opt, 0}), ptr_to_circle_event({nullptr, nullptr, h2})));                    
                }
                else if(p.x < q2.x){
                    std::optional<point> q2_opt = q2;
                    if(q1_opt.has_value()){
                        h1.start_point_x = (p.x + q1_opt.value().x)/2; h1.start_point_y = std::nullopt; h1.u = 0; h1.v = -1;
                    }
                    h2.start_point_x = (p.x + q2.x)/2; h2.start_point_y = std::nullopt; h2.u = 0; h2.v = -1;
                    beachline.erase(it1);
                    beachline.insert(std::make_pair(node({q1_opt, p_opt, 0}), ptr_to_circle_event({nullptr, nullptr, h1})));
                    beachline.insert(std::make_pair(node({p_opt, q2_opt, 0}), ptr_to_circle_event({nullptr, nullptr, h2})));
                }
                else{
                    // std::cout<< "p.x = " << p.x << " q2.x = " << q2.x << " q3.x = " << q3.x << std::endl;
                    // std::cout<< "ERROR: Site event at degenerate stage\n" << std::endl;
                    // assert(false);
                }
            }
        }
        return;
    }
    
    // Find the arc above the new point
    auto it1 = beachline.lower_bound(node({p, std::nullopt, 1}));
    auto it2 = beachline.upper_bound(node({p, std::nullopt, 1}));

    // Since *beachline.begin() will be of the form (NULL, q) and *(--beachline.end()) will be of the form (p, NULL),
    // the new point will always be between these two points
    // assert (it1 != beachline.begin() && it1 != beachline.end());
    // assert (it2 != beachline.begin() && it2 != beachline.end());

    if(it1 == it2){     // The new point is between two breakpoints

        // Checking if the arc has a pointer to a circle event in Q, then this circle event is a false alarm.
        // The false circle event must be deleted from Q.
        if(it1->second.p_circle_ptr){
            auto it3 = it1;
            it3--;
            // assert(it1->second.p_circle_ptr->point_type == 1 && it3->second.q_circle_ptr->point_type == 1);
            event_queue.erase(event_queue.find(*(it1->second.p_circle_ptr)));
            free(it1->second.p_circle_ptr);
            it1->second.p_circle_ptr = nullptr;
            it3->second.q_circle_ptr = nullptr;
        }
        // assert(it1->first.p_opt.has_value());
        point q1 = it1->first.p_opt.value();
        it1--;
        // assert(it1->first.q_opt.has_value());
        point q2 = it1->first.q_opt.value();
        // assert(q1.x == q2.x and q1.y == q2.y);  // (q1 = q2) At this stage, the nodes are in order it1-(q3, q), it2-(q, q4) where (it1, it2) is the arc above p

        std::optional<point> q_opt = q1;
        std::optional<point> p_opt = p;

        // Check the triple of consecutive arcs where the new arc for point p is the left arc to see if the breakpoints converge. 
        // If so, insert the circle event into event queue and add pointers between the node in beachline and the node in event queue. 
        // Do the same for the triple where the new arc is the right arc.
        
        std::pair<half_edge, half_edge> half_edge_pair = start_edge(q1, p.x, sweepLine_y);
        half_edge h_left = half_edge_pair.first;
        half_edge h_right = half_edge_pair.second;

        node n_left = node({q_opt, p_opt, 0});
        auto it_left = beachline.insert(std::make_pair(n_left, ptr_to_circle_event({nullptr, nullptr, h_left})));
        if(it1->first.p_opt.has_value()){
            point q3 = it1->first.p_opt.value();
            node n_left_left = it1->first;
            std::optional<point> left_circle_point = findCircleEvent(q3, q2, p);
            if(left_circle_point.has_value()){
                node* n_left_ptr = new node;
                *n_left_ptr = n_left;
                node* n_left_left_ptr = new node;
                *n_left_left_ptr = n_left_left;
                point* l_point_ptr = new point; 
                point l_point = left_circle_point.value();
                *l_point_ptr = l_point;
                l_point_ptr->point_type = 1;
                it_left->second.p_circle_ptr = l_point_ptr;
                it1->second.q_circle_ptr = l_point_ptr;
                event_queue.insert(std::make_pair(l_point, ptr_to_arc({n_left_left_ptr, n_left_ptr})));
            }
            else{
                it_left->second.p_circle_ptr = nullptr;
                it1->second.q_circle_ptr = nullptr;
            }
        }

        node n_right = node({p_opt, q_opt, 0});
        auto it_right = beachline.insert(std::make_pair(n_right, ptr_to_circle_event({nullptr, nullptr, h_right})));
        if(it2->first.q_opt.has_value()){
            point q4 = it2->first.q_opt.value();
            node n_right_right = it2->first;
            std::optional<point> right_circle_point = findCircleEvent(p, q1, q4);
            if(right_circle_point.has_value()){
                node* n_right_ptr = new node;
                *n_right_ptr = n_right;
                node* n_right_right_ptr = new node;
                *n_right_right_ptr = n_right_right;
                point* r_point_ptr = new point;
                point r_point = right_circle_point.value();
                *r_point_ptr = r_point;
                r_point_ptr->point_type = 1;
                it_right->second.q_circle_ptr = r_point_ptr;
                it2->second.p_circle_ptr = r_point_ptr;
                event_queue.insert(std::make_pair(r_point, ptr_to_arc({n_right_ptr, n_right_right_ptr})));
            }
            else{
                it_right->second.q_circle_ptr = nullptr;
                it2->second.p_circle_ptr = nullptr;
            }
        }        
    }
    else{       // The new point is directly below a breakpoint
        // std::cout<< "ERROR: Site event at degenerate stage\n" << std::endl;
        // assert(it1->first.p_opt.has_value() && it1->first.q_opt.has_value() && it2->first.p_opt.has_value());
        point q1 = it1->first.p_opt.value();
        point q2 = it1->first.q_opt.value();
        point q3 = it2->first.p_opt.value();
        // assert(q2.x == q3.x and q2.y == q3.y);

        half_edge h_old = it1->second.half_edge_traced;
        std::optional<point> circumcentre_opt = findCircleEvent(p,q1,q2);
        // assert(circumcentre_opt.has_value());
        point circumcentre = circumcentre_opt.value();
        h_old.end_point_x = circumcentre.x; h_old.end_point_y = circumcentre.y;
        half_edges.push_back(h_old);

        std::optional<point> q1_opt = q1;
        std::optional<point> q2_opt = q2;
        std::optional<point> p_opt = p;
    
        // it1 = breakpoint directly above p
        // it2 = breakpoint to the right of p

        // If adjacent arcs have circle events in Q, then remove them. (False alarms)
        // Left arc - Checking for left circle event and removing it from event queue if it exists
        if(it1->second.p_circle_ptr){
            auto it3 = it1;
            it3--;
            // assert(it1->second.p_circle_ptr->point_type == 1 && it3->second.q_circle_ptr->point_type == 1);
            event_queue.erase(event_queue.find(*(it1->second.p_circle_ptr)));
            free(it1->second.p_circle_ptr);
            it1->second.p_circle_ptr = nullptr;
            it3->second.q_circle_ptr = nullptr;
        }

        // Right arc - Checking for right circle event and removing it from event queue if it exists
        if(it1->second.q_circle_ptr){
            auto it3 = it1;
            it3++;
            // assert(it1->second.q_circle_ptr->point_type == 1 && it3->second.p_circle_ptr->point_type == 1);
            event_queue.erase(event_queue.find(*(it1->second.q_circle_ptr)));
            free(it1->second.p_circle_ptr);
            it1->second.q_circle_ptr = nullptr;
            it3->second.p_circle_ptr = nullptr;
        }

        // Erases it1 (breakpoint above p) and inserts the two new breakpoints
        beachline.erase(it1);
        half_edge h_left = start_edge(q1, p.x, sweepLine_y).first;
        half_edge h_right = start_edge(q2, p.x, sweepLine_y).second;
        node n_left = node({q1_opt, p_opt, 0});
        auto it_left = beachline.insert(std::make_pair(n_left, ptr_to_circle_event({nullptr, nullptr, h_left})));
        node n_right = node({p_opt, q2_opt, 0});
        auto it_right = beachline.insert(std::make_pair(n_right, ptr_to_circle_event({nullptr, nullptr, h_right})));

        // Soon, the nodes will be in order it_left-(q1, p), it_right-(p, q2) where (it_left, it_right) is the arc inserted

        // Check the triple of consecutive arcs where the new arc for point p is the left arc to see if the breakpoints converge. 
        // If so, insert the circle event into event queue and add pointers between the node in beachline and the node in event queue. 
        // Do the same for the triple where the new arc is the right arc.
        
        auto it_left_left = it_left; it_left_left--;
        if(it_left_left->first.p_opt.has_value()){
            point p_left_left = it_left_left->first.p_opt.value();
            node n_left_left = it_left_left->first;
            std::optional<point> left_circle_point = findCircleEvent(p_left_left, q1, p);
            if(left_circle_point.has_value()){
                node* n_left_ptr = new node;
                *n_left_ptr = n_left;
                node* n_left_left_ptr = new node;
                *n_left_left_ptr = n_left_left;
                point* l_point_ptr = new point; 
                point l_point = left_circle_point.value();
                *l_point_ptr = l_point;
                l_point_ptr->point_type = 1;
                it_left->second.p_circle_ptr = l_point_ptr;
                it_left_left->second.q_circle_ptr = l_point_ptr;
                event_queue.insert(std::make_pair(l_point, ptr_to_arc({n_left_left_ptr, n_left_ptr})));
            }
            else{
                it_left->second.p_circle_ptr = nullptr;
                it_left_left->second.q_circle_ptr = nullptr;
            }
        }

        auto it_right_right = it_right; it_right_right++;
        if(it_right_right->first.q_opt.has_value()){
            point q_right_right = it_right_right->first.q_opt.value();
            node n_right_right = it_right_right->first;
            std::optional<point> right_circle_point = findCircleEvent(p, q2, q_right_right);
            if(right_circle_point.has_value()){
                node* n_right_ptr = new node;
                *n_right_ptr = n_right;
                node* n_right_right_ptr = new node;
                *n_right_right_ptr = n_right_right;
                point* r_point_ptr = new point;
                point r_point = right_circle_point.value();
                *r_point_ptr = r_point;
                r_point_ptr->point_type = 1;
                it_right->second.q_circle_ptr = r_point_ptr;
                it_right_right->second.p_circle_ptr = r_point_ptr;
                event_queue.insert(std::make_pair(r_point, ptr_to_arc({n_right_ptr, n_right_right_ptr})));
            }
            else{
                it_right->second.q_circle_ptr = nullptr;
                it_right_right->second.p_circle_ptr = nullptr;
            }
        }
    }
}

// Function to handle the circle event
void HandleCircleEvent(node left_node, node right_node){
    // assert(left_node.isEventPoint == 0 && right_node.isEventPoint == 0);
    // assert(left_node.p_opt.has_value() && left_node.q_opt.has_value() && right_node.p_opt.has_value() && right_node.q_opt.has_value());
    
    // std::cout << "\nInside HandleCircleEvent(), " << std::endl;
    // std::cout<< "Left node: " << left_node.p_opt.has_value() << " " << left_node.q_opt.has_value() << " " << left_node.isEventPoint << std::endl;
    // if(left_node.p_opt.has_value()) std::cout<< "Left node: " << left_node.p_opt.value().index << " " << left_node.p_opt.value().x << " " << left_node.p_opt.value().y << std::endl;
    // if(left_node.q_opt.has_value()) std::cout<< "Left node: " << left_node.q_opt.value().index << " " << left_node.q_opt.value().x << " " << left_node.q_opt.value().y << std::endl;
    // std::cout<< "Right node: " << right_node.p_opt.has_value() << " " << right_node.q_opt.has_value() << " " << right_node.isEventPoint << std::endl;
    // if(right_node.p_opt.has_value()) std::cout<< "Right node: " << right_node.p_opt.value().index << " " << right_node.p_opt.value().x << " " << right_node.p_opt.value().y << std::endl;
    // if(right_node.q_opt.has_value()) std::cout<< "Right node: " << right_node.q_opt.value().index << " " << right_node.q_opt.value().x << " " << right_node.q_opt.value().y << std::endl;

    sweepLine_y += epsilon2;
    std::cout << "Sweepline : Y = " << sweepLine_y << std::endl;
    auto it2 = find_custom(left_node);
    auto it3 = find_custom(right_node);
    // assert(it2 != beachline.begin()); // assert(it3 != beachline.begin());
    // assert(it2 != beachline.end()); // assert(it3 != beachline.end());
    // assert(it2->first.q_opt.has_value() && it3->first.p_opt.has_value());

    // std::cout<< "(" << it2->first.q_opt.value().x <<"," << it2->first.q_opt.value().y << ") (" << it3->first.p_opt.value().x <<"," << it3->first.p_opt.value().y <<")" << std::endl;
    sweepLine_y -= epsilon2;
    // assert(it2->second.q_circle_ptr && it3->second.p_circle_ptr);
    // assert(it2->second.q_circle_ptr->point_type==1 && it3->second.p_circle_ptr->point_type==1);
    auto it1 = it2; it1--;
    auto it4 = it3; it4++;

    // Now, the iterators point to the breakpoints in the order (it1-(p1,p2), it2-(p2,p3), it3-(p3,p4), it4-(p4,p5)) where (it2, it3) is the arc to be removed
    // assert(it1->first.q_opt.has_value() && it2->first.p_opt.has_value());
    // assert(it1->first.q_opt.value().x == it2->first.p_opt.value().x && it1->first.q_opt.value().y == it2->first.p_opt.value().y);
    point p2 = it2->first.p_opt.value();
    // assert(it2->first.q_opt.has_value() && it3->first.p_opt.has_value());
    // assert(it2->first.q_opt.value().x == it3->first.p_opt.value().x && it2->first.q_opt.value().y == it3->first.p_opt.value().y);
    point p3 = it3->first.p_opt.value();
    // assert(it3->first.q_opt.has_value() && it4->first.p_opt.has_value());
    // std::cout<< it3->first.q_opt.value().x << " " << it3->first.q_opt.value().y << std::endl;
    // std::cout<< it4->first.p_opt.value().x << " " << it4->first.p_opt.value().y << std::endl;
    // assert(it3->first.q_opt.value().x == it4->first.p_opt.value().x && it3->first.q_opt.value().y == it4->first.p_opt.value().y);
    point p4 = it4->first.p_opt.value();

    // Checking for left circle event and removing it from event queue if it exists
    if(it2->second.p_circle_ptr){
        // assert(it2->second.p_circle_ptr->point_type == 1 && it1->second.q_circle_ptr->point_type == 1);
        event_queue.erase(event_queue.find(*(it2->second.p_circle_ptr)));
        free(it2->second.p_circle_ptr);
        it1->second.q_circle_ptr = nullptr;
        it2->second.p_circle_ptr = nullptr;
    }

    // Checking for right circle event and removing it from event queue if it exists
    if(it3->second.q_circle_ptr){
        // assert(it3->second.q_circle_ptr->point_type == 1 && it4->second.p_circle_ptr->point_type == 1);
        event_queue.erase(event_queue.find(*(it3->second.q_circle_ptr)));
        free(it3->second.q_circle_ptr);
        it3->second.q_circle_ptr = nullptr;
        it4->second.p_circle_ptr = nullptr;
    }

    std::optional<point> circumcentre_opt = findCircumcentre(p2, p3, p4);
    // assert(circumcentre_opt.has_value());
    point circumcentre = circumcentre_opt.value();
    half_edge h_left = it2->second.half_edge_traced;
    half_edge h_right = it3->second.half_edge_traced;
    h_left.end_point_x = circumcentre.x; h_left.end_point_y = circumcentre.y;
    h_right.end_point_x = circumcentre.x; h_right.end_point_y = circumcentre.y;
    half_edges.push_back(h_left);
    half_edges.push_back(h_right);
    // Erase two middle breakpoints and insert the new breakpoint into the beachline
    beachline.erase(it2);
    beachline.erase(it3);
    // std::cout<< "After erasing the two merging breakpoints : \n";
    // printBeachline();

    node* n_ptr = new node;
    node n = node({p2, p4, 0});
    *n_ptr = n;
    half_edge h_new = half_edge();
    h_new.start_point_x = circumcentre.x; h_new.start_point_y = circumcentre.y; h_new.u = (p4.y - p2.y); h_new.v = (p2.x - p4.x);

    auto it = beachline.insert(std::make_pair(n, ptr_to_circle_event({nullptr, nullptr, h_new})));
    // std::cout<< "After insertion of new breakpoint : \n";
    // printBeachline();

    // Check triples of new consecutive events for new circle events and insert them into the event queue
    if(it1->first.p_opt.has_value()){
        point p1 = it1->first.p_opt.value();
        std::optional<point> left_circle_point = findCircleEvent(p1, p2, p4);
        if(left_circle_point.has_value()){
            // std::cout<< "Left circle event: " << left_circle_point.value().x << " " << left_circle_point.value().y << std::endl;
            node n1 = it1->first;
            node* n1_ptr = new node;
            *n1_ptr = n1;
            point l_point = left_circle_point.value();
            point* l_point_ptr = new point; 
            *l_point_ptr = l_point;
            l_point_ptr->point_type = 1;
            it->second.p_circle_ptr = l_point_ptr;
            it1->second.q_circle_ptr = l_point_ptr;
            event_queue.insert(std::make_pair(*l_point_ptr, ptr_to_arc({n1_ptr, n_ptr})));
        }
        else{
            // std::cout<< "No left circle event" << std::endl;
            it->second.p_circle_ptr = nullptr;
            it1->second.q_circle_ptr = nullptr;
        }
    }

    if(it4->first.q_opt.has_value()){
        point p5 = it4->first.q_opt.value();
        std::optional<point> right_circle_point = findCircleEvent(p2, p4, p5);
        if(right_circle_point.has_value()){
            // std::cout<< "Right circle event: " << right_circle_point.value().x << " " << right_circle_point.value().y << std::endl;
            node n4 = it4->first;
            node* n4_ptr = new node;
            *n4_ptr = n4;
            point r_point = right_circle_point.value();
            point* r_point_ptr = new point;
            *r_point_ptr = r_point;
            r_point_ptr->point_type = 1;
            it->second.q_circle_ptr = r_point_ptr;
            it4->second.p_circle_ptr = r_point_ptr;
            event_queue.insert(std::make_pair(*r_point_ptr, ptr_to_arc({n_ptr, n4_ptr})));
        }
        else{
            // std::cout<< "No right circle event" << std::endl;
            it->second.q_circle_ptr = nullptr;
            it4->second.p_circle_ptr = nullptr;
        }
    }
}

int main(){
    int n; std::cout<< "Number of points: "; std::cin >> n;
    double x,y;
    std::cout<< "Enter the points: \n";

	for(int i=0; i<n; i++){
        std::cin >> x >> y;
        input_points.push_back(point({i, x, y, 0}));
		event_queue.insert(std::make_pair(point({i, x, y, 0}), ptr_to_arc()));	
	}

    while(!event_queue.empty()){
        point p = event_queue.begin()->first;
        ptr_to_arc arc = event_queue.begin()->second;

        event_queue.erase(event_queue.begin());
        sweepLine_y = p.y;

        std::cout<< "\n---------------------------------------------" << std::endl;
        std::cout << "\nPROCESSING POINT: " << p.index << " (" << p.x << ", " << p.y << ") " << p.point_type << '\n' << std::endl;

        if(p.point_type == 0){  // Site event
            HandleSiteEvent(p);
        }
        else{   // Circle event
            // printBeachline();
            // assert(arc.left_end && arc.right_end);
            node left_node = *arc.left_end;
            node right_node = *arc.right_end;
            HandleCircleEvent(left_node, right_node);
        }

        // Debugging purposes
        // printBeachline();
        // printEventQueue();
        // printHalfEdges();
        // std::cout<< "\n---------------------------------------------" << std::endl;
    }
    
    for(auto it = beachline.begin(); it != beachline.end(); it++){
        if(!(it->second.half_edge_traced.u == 0 && it->second.half_edge_traced.v == 0)){
            // assert(!it->second.half_edge_traced.end_point_x.has_value() && !it->second.half_edge_traced.end_point_y.has_value());
            half_edges.push_back(it->second.half_edge_traced);
        }
    }
    // printHalfEdges();
    writeToFile(half_edges);
    return 0;
}