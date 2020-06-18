//
//  main.cpp
//  User-centric NOMA
//
//  Created by 俞斩味 on 2020/3/20.
//  Copyright © 2020 俞斩味. All rights reserved.
//

#define BS 21
#define UE 210
#define bandwidth 180000
#define noise 7.16593e-16
#define sP 50e-3
#define mP 200e-3
#define demand 0.03
#define rho_limit 1

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>

using namespace std;

double gain[BS][UE];
double association[BS][UE];
double RU[BS];
double Power[BS];
double alpha[UE];
int is_selected[UE] = {0};//a flag showing UE if is selected
int ndx;//send i to global

void initialize()
{
    for (int i = 0; i < BS; ++i)
        RU[i] = 0;
    for (int j = 0; j < UE; ++j)
        alpha[j] = 1;
    for (int j = 0; j < UE; ++j)
        is_selected[j] = 0;
}

void read()//read matrix
{
    for (int i = 0; i < 21; ++i)
        RU[i] = 1;
    ifstream gain_matrix("/Users/ZwYu/Desktop/Codes/User-centric NOMA/gain_matrix");
    ifstream association_matrix("/Users/ZwYu/Desktop/Codes/User-centric NOMA/cell_user_matrix");
    string buff;
    int k = 0;
    while (getline(gain_matrix,buff,','))
    {
        if ((k+1)%210 == 0)
        {
            *(gain[k/210]+(k%210)) = atof((buff.substr(0,17)).c_str());
            k++;
            *(gain[k/210]+(k%210)) = atof((buff.substr(18,17)).c_str());
        }
        else
        {
            *(gain[k/210]+(k%210)) = atof(buff.c_str());
        }
        k++;
    }
    string buff2;
    int p = 0;
    while (getline(association_matrix,buff2,','))
    {
        if ((p+1)%210 == 0)
        {
            *(association[p/210]+(p%210)) = atof((buff2.substr(0,1)).c_str());
            p++;
            *(association[p/210]+(p%210)) = atof((buff2.substr(2,1)).c_str());
        }
        else
        {
            *(association[p/210]+(p%210)) = atof(buff2.c_str());
        }
        p++;
    }
}

void converted()//convert scaled
{
    for (int i = 0; i < 21; ++i)
        for (int j = 0; j < 210; ++j)
            gain[i][j] = pow(10,gain[i][j]/10);
    for (int i = 0; i < 7; ++i)
        Power[i] = mP;
    for (int i = 7; i < 21; ++i)
        Power[i] = sP;
}

double calculate_xu(int j, int h, int i)//j and h is UE, i is both associated BS
{
    double w1,w2;
    int worse;//to record the worse UE index
    
    double wj = 0;
    for (int k = 0; k < BS; ++k)
        wj += (!association[k][j])*Power[k]*gain[k][j]*RU[k];
    wj = (wj + noise)/gain[i][j];
    
    double wh = 0;
    for (int k = 0; k < BS; ++k)
        wh += (!association[k][h])*Power[k]*gain[k][h]*RU[k];
    wh = (wh + noise)/gain[i][h];
    
    if (wj >= wh)
    {
        w1 = wh;
        w2 = wj;
        worse = j;
    }
    else
    {
        w1 = wj;
        w2 = wh;
        worse = h;
    }
    
    double xu_min = 0.0000001;
    double xu_max = 10;
    double xu_current = 0;
    while ((xu_max - xu_min) > 0.00000001)
    {
        xu_current = (xu_max + xu_min) / 2;
        double P = w1 * pow(2, (alpha[j]+alpha[h])*demand/xu_current) + (w2 - w1) * pow(2, alpha[worse]*demand/xu_current) - w2;
        if (P > Power[i])
            xu_min = xu_current;
        else
            xu_max = xu_current;
    }
    return xu_current;
}

double calculate_x(int j, int i)
{
    double IN = 0;
    for (int k = 0; k < BS; ++k)
    {
        if (k != i)
            IN += Power[k] * gain[k][j] * RU[k];
    }
    IN += noise;
    double x = alpha[j]*demand/log(1+Power[i]*gain[i][j]/IN)*log(2);
    return x;
}


typedef long long s64;

const int INF = 2147483647;

const int MaxN = 22;
const int MaxM = 79800;

template <class T>
inline void tension(T &a, const T &b)
{
    if (b < a)
        a = b;
}
template <class T>
inline void relax(T &a, const T &b)
{
    if (b > a)
        a = b;
}
template <class T>
inline int size(const T &a)
{
    return (int)a.size();
}

inline int getint()
{
    char c;
    while (c = getchar(), '0' > c || c > '9');

    int res = c - '0';
    while (c = getchar(), '0' <= c && c <= '9')
        res = res * 10 + c - '0';
    return res;
}

const int MaxNX = MaxN + MaxN;

struct edge
{
    int v, u, w;

    edge(){}
    edge(const int &_v, const int &_u, const int &_w)
        : v(_v), u(_u), w(_w){}
};

int n, m;
edge mat[MaxNX + 1][MaxNX + 1];

int n_matches;
s64 tot_weight;
int mate[MaxNX + 1];
int lab[MaxNX + 1];

int q_n, q[MaxN];
int fa[MaxNX + 1], col[MaxNX + 1];
int slackv[MaxNX + 1];

int n_x;
int bel[MaxNX + 1], blofrom[MaxNX + 1][MaxN + 1];
vector<int> bloch[MaxNX + 1];

inline int e_delta(const edge &e) // does not work inside blossoms
{
    return lab[e.v] + lab[e.u] - mat[e.v][e.u].w * 2;
}
inline void update_slackv(int v, int x)
{
    if (!slackv[x] || e_delta(mat[v][x]) < e_delta(mat[slackv[x]][x]))
        slackv[x] = v;
}
inline void calc_slackv(int x)
{
    slackv[x] = 0;
    for (int v = 1; v <= n; v++)
        if (mat[v][x].w > 0 && bel[v] != x && col[bel[v]] == 0)
            update_slackv(v, x);
}

inline void q_push(int x)
{
    if (x <= n)
        q[q_n++] = x;
    else
    {
        for (int i = 0; i < size(bloch[x]); i++)
            q_push(bloch[x][i]);
    }
}
inline void set_mate(int xv, int xu)
{
    mate[xv] = mat[xv][xu].u;
    if (xv > n)
    {
        edge e = mat[xv][xu];
        int xr = blofrom[xv][e.v];
        int pr = find(bloch[xv].begin(), bloch[xv].end(), xr) - bloch[xv].begin();
        if (pr % 2 == 1)
        {
            reverse(bloch[xv].begin() + 1, bloch[xv].end());
            pr = size(bloch[xv]) - pr;
        }

        for (int i = 0; i < pr; i++)
            set_mate(bloch[xv][i], bloch[xv][i ^ 1]);
        set_mate(xr, xu);

        rotate(bloch[xv].begin(), bloch[xv].begin() + pr, bloch[xv].end());
    }
}
inline void set_bel(int x, int b)
{
    bel[x] = b;
    if (x > n)
    {
        for (int i = 0; i < size(bloch[x]); i++)
            set_bel(bloch[x][i], b);
    }
}

inline void augment(int xv, int xu)
{
    while (true)
    {
        int xnu = bel[mate[xv]];
        set_mate(xv, xu);
        if (!xnu)
            return;
        set_mate(xnu, bel[fa[xnu]]);
        xv = bel[fa[xnu]], xu = xnu;
    }
}
inline int get_lca(int xv, int xu)
{
    static bool book[MaxNX + 1];
    for (int x = 1; x <= n_x; x++)
        book[x] = false;
    while (xv || xu)
    {
        if (xv)
        {
            if (book[xv])
                return xv;
            book[xv] = true;
            xv = bel[mate[xv]];
            if (xv)
                xv = bel[fa[xv]];
        }
        swap(xv, xu);
    }
    return 0;
}

inline void add_blossom(int xv, int xa, int xu)
{
    int b = n + 1;
    while (b <= n_x && bel[b])
        b++;
    if (b > n_x)
        n_x++;

    lab[b] = 0;
    col[b] = 0;

    mate[b] = mate[xa];

    bloch[b].clear();
    bloch[b].push_back(xa);
    for (int x = xv; x != xa; x = bel[fa[bel[mate[x]]]])
        bloch[b].push_back(x), bloch[b].push_back(bel[mate[x]]), q_push(bel[mate[x]]);
    reverse(bloch[b].begin() + 1, bloch[b].end());
    for (int x = xu; x != xa; x = bel[fa[bel[mate[x]]]])
        bloch[b].push_back(x), bloch[b].push_back(bel[mate[x]]), q_push(bel[mate[x]]);

    set_bel(b, b);

    for (int x = 1; x <= n_x; x++)
    {
        mat[b][x].w = mat[x][b].w = 0;
        blofrom[b][x] = 0;
    }
    for (int i = 0; i < size(bloch[b]); i++)
    {
        int xs = bloch[b][i];
        for (int x = 1; x <= n_x; x++)
            if (mat[b][x].w == 0 || e_delta(mat[xs][x]) < e_delta(mat[b][x]))
                mat[b][x] = mat[xs][x], mat[x][b] = mat[x][xs];
        for (int x = 1; x <= n_x; x++)
            if (blofrom[xs][x])
                blofrom[b][x] = xs;
    }
    calc_slackv(b);
}
inline void expand_blossom1(int b) // lab[b] == 1
{
    for (int i = 0; i < size(bloch[b]); i++)
        set_bel(bloch[b][i], bloch[b][i]);

    int xr = blofrom[b][mat[b][fa[b]].v];
    int pr = find(bloch[b].begin(), bloch[b].end(), xr) - bloch[b].begin();
    if (pr % 2 == 1)
    {
        reverse(bloch[b].begin() + 1, bloch[b].end());
        pr = size(bloch[b]) - pr;
    }

    for (int i = 0; i < pr; i += 2)
    {
        int xs = bloch[b][i], xns = bloch[b][i + 1];
        fa[xs] = mat[xns][xs].v;
        col[xs] = 1, col[xns] = 0;
        slackv[xs] = 0, calc_slackv(xns);
        q_push(xns);
    }
    col[xr] = 1;
    fa[xr] = fa[b];
    for (int i = pr + 1; i < size(bloch[b]); i++)
    {
        int xs = bloch[b][i];
        col[xs] = -1;
        calc_slackv(xs);
    }

    bel[b] = 0;
}
inline void expand_blossom_final(int b) // at the final stage
{
    for (int i = 0; i < size(bloch[b]); i++)
    {
        if (bloch[b][i] > n && lab[bloch[b][i]] == 0)
            expand_blossom_final(bloch[b][i]);
        else
            set_bel(bloch[b][i], bloch[b][i]);
    }
    bel[b] = 0;
}

inline bool on_found_edge(const edge &e)
{
    int xv = bel[e.v], xu = bel[e.u];
    if (col[xu] == -1)
    {
        int nv = bel[mate[xu]];
        fa[xu] = e.v;
        col[xu] = 1, col[nv] = 0;
        slackv[xu] = slackv[nv] = 0;
        q_push(nv);
    }
    else if (col[xu] == 0)
    {
        int xa = get_lca(xv, xu);
        if (!xa)
        {
            augment(xv, xu), augment(xu, xv);
            for (int b = n + 1; b <= n_x; b++)
                if (bel[b] == b && lab[b] == 0)
                    expand_blossom_final(b);
            return true;
        }
        else
            add_blossom(xv, xa, xu);
    }
    return false;
}

bool match()
{
    for (int x = 1; x <= n_x; x++)
        col[x] = -1, slackv[x] = 0;

    q_n = 0;
    for (int x = 1; x <= n_x; x++)
        if (bel[x] == x && !mate[x])
            fa[x] = 0, col[x] = 0, slackv[x] = 0, q_push(x);
    if (q_n == 0)
        return false;

    while (true)
    {
        for (int i = 0; i < q_n; i++)
        {
            int v = q[i];
            for (int u = 1; u <= n; u++)
                if (mat[v][u].w > 0 && bel[v] != bel[u])
                {
                    int d = e_delta(mat[v][u]);
                    if (d == 0)
                    {
                        if (on_found_edge(mat[v][u]))
                            return true;
                    }
                    else if (col[bel[u]] == -1 || col[bel[u]] == 0)
                        update_slackv(v, bel[u]);
                }
        }

        int d = INF;
        for (int v = 1; v <= n; v++)
            if (col[bel[v]] == 0)
                tension(d, lab[v]);
        for (int b = n + 1; b <= n_x; b++)
            if (bel[b] == b && col[b] == 1)
                tension(d, lab[b] / 2);
        for (int x = 1; x <= n_x; x++)
            if (bel[x] == x && slackv[x])
            {
                if (col[x] == -1)
                    tension(d, e_delta(mat[slackv[x]][x]));
                else if (col[x] == 0)
                    tension(d, e_delta(mat[slackv[x]][x]) / 2);
            }

        for (int v = 1; v <= n; v++)
        {
            if (col[bel[v]] == 0)
                lab[v] -= d;
            else if (col[bel[v]] == 1)
                lab[v] += d;
        }
        for (int b = n + 1; b <= n_x; b++)
            if (bel[b] == b)
            {
                if (col[bel[b]] == 0)
                    lab[b] += d * 2;
                else if (col[bel[b]] == 1)
                    lab[b] -= d * 2;
            }

        q_n = 0;
        for (int v = 1; v <= n; v++)
            if (lab[v] == 0) // all unmatched vertices' labels are zero! cheers!
                return false;
        for (int x = 1; x <= n_x; x++)
            if (bel[x] == x && slackv[x] && bel[slackv[x]] != x && e_delta(mat[slackv[x]][x]) == 0)
            {
                if (on_found_edge(mat[slackv[x]][x]))
                    return true;
            }
        for (int b = n + 1; b <= n_x; b++)
            if (bel[b] == b && col[b] == 1 && lab[b] == 0)
                expand_blossom1(b);
    }
    return false;
}

void calc_max_weight_match()
{
    for (int v = 1; v <= n; v++)
        mate[v] = 0;

    n_x = n;
    n_matches = 0;
    tot_weight = 0;

    bel[0] = 0;
    for (int v = 1; v <= n; v++)
        bel[v] = v, bloch[v].clear();
    for (int v = 1; v <= n; v++)
        for (int u = 1; u <= n; u++)
            blofrom[v][u] = v == u ? v : 0;

    int w_max = 0;
    for (int v = 1; v <= n; v++)
        for (int u = 1; u <= n; u++)
            relax(w_max, mat[v][u].w);
    for (int v = 1; v <= n; v++)
        lab[v] = w_max;

    while (match())
        n_matches++;

    for (int v = 1; v <= n; v++)
        if (mate[v] && mate[v] < v)
            tot_weight += mat[v][mate[v]].w;
}

/*int main()//MWM
{
    cin>> n>> m;

    for (int v = 1; v <= n; v++)
        for (int u = 1; u <= n; u++)
            mat[v][u] = edge(v, u, 0);

    for (int i = 0; i < m; i++)
    {
        int v, u, w;
        cin>> v>> u>> w;
        mat[v][u].w = mat[u][v].w = w;
    }

    calc_max_weight_match();

    cout<<tot_weight<<"\n";
    for (int v = 1; v <= n; v++)
        cout<<mate[v]<<" ";
    cout<<"\n";

    return 0;
}*/

double MWM (int i)
{
    int m_count = 0;
    int n_count = 0;
    for (int j = 0; j < UE; ++j)
    {
        if (association[i][j])
            {
                n_count++;//tongji n
                for (int h = 0; h < UE; ++h)
                    if (association[i][h] && (j!=h) && (j<h))
                    m_count++;//tongji m
            }
    }
    //cout<<n_count<<"\n";
    //cout<<m_count<<"\n";
    
    if (n_count % 2 == 0)
    {
        n = n_count;
        m = m_count;
    }
    else
    {
        n = n_count + 1;
        m = m_count + n;
    }

    for (int v = 1; v <= n; v++)//clear cache
        for (int u = 1; u <= n; u++)
            mat[v][u] = edge(v, u, 0);

    int v, u;
    v = 1;
    for (int j = 0; j < UE; ++j)//read graph
    {
        u = v+1;
        if (association[i][j])
            {
                if (n_count % 2 == 1)
                {
                    mat[n][v].w = mat[v][n].w = (int)((10-calculate_x(j, i))*1000000);
                }
                for (int h = 0; h < UE; ++h)
                    if (association[i][h] && (j<h))
                    {
                        mat[v][u].w = mat[u][v].w = (int)((10-calculate_xu(j, h, i))*1000000);
                        //cout<<v<<" "<<u<<" "<<mat[v][u].w<<"\n";
                        u++;
                    }
                v++;
            }
    }

    calc_max_weight_match();
    
    /*cout<<tot_weight<<"\n";
    for (int v = 1; v <= n; v++)
        cout<<mate[v]<<" ";
    cout<<"\n";*/
    
    return 10*n/2-tot_weight/1e6;//sigma xu

}

double calculate_sum()
{
    double sum = 0;
    for (int i = 0; i < BS; ++i)
        sum += RU[i];
    return sum;
}

bool cmp(int a, int b)
{
    return gain[ndx][a] < gain[ndx][b];
}

void select_node(double percent, double given_alpha)
{
    int BS_count_UE[BS] = {0};//count how much UEs are associted to BSs
    int BS_selected_number[BS];//show number of selected UE
    int current_BS_association[UE] = {0};
    for (int i = 0; i < BS; ++i)
        for (int j = 0; j < UE; ++j)
            if (association[i][j])
                BS_count_UE[i]++;
    for (int i = 0; i < BS; ++i)
        BS_selected_number[i] = (int) (BS_count_UE[i] * percent);
    
    for (int i = 0; i < BS; ++i)
    {
        int n = 0;
        for (int j = 0; j < UE; ++j)
            if (association[i][j])
                current_BS_association[n++] = j;
        ndx = i;
        sort(current_BS_association, current_BS_association + BS_count_UE[i], cmp);
        for (int j = 0; j < BS_selected_number[i]; ++j)
            is_selected[current_BS_association[j]] = 1;
    }
    for (int j = 0; j < UE; ++j)
        if (is_selected[j])
            alpha[j] = given_alpha;
        else
            alpha[j] = 1;
}

double get_maximum_RU(double given_percent, double given_alpha)
{
    initialize();
    
    select_node(given_percent, given_alpha);
    
    double RU_k[BS];
    
    int condition_no_converge_flag = 1;
    while (condition_no_converge_flag)
    {
        condition_no_converge_flag = 0;
        for (int i = 0; i < BS; ++i)
            RU_k[i] = RU[i];
        for (int i = 0; i < BS; ++i)
            RU[i] = MWM(i);
        for (int i = 0; i < BS; ++i)
            if (abs(RU_k[i] - RU[i]) > 0.0001)
                condition_no_converge_flag++;
    }
    
    sort(RU, RU+BS);
    for (int i = 0; i < BS; ++i)
    {
        //cout << i+1 << ":" << RU[i] << "\n";
    }
    //cout<<"Total:"<<calculate_sum()<<"\n";
    return RU[BS-1];
}


int main()
{
    read();
    association[1][43] = 0;
    converted();
    
    double alpha_max = 10;
    double alpha_min = 1;
    double alpha_avg;
    double maximum_RU;
    int count = 0;
    while (alpha_max - alpha_min > 0.000001)
    {
        alpha_avg = (alpha_max + alpha_min) / 2;
        maximum_RU = get_maximum_RU(1, alpha_avg);
        if (maximum_RU > rho_limit)
            alpha_max = alpha_avg;
        else
            alpha_min = alpha_avg;
        count++;
    }
    cout<<alpha_avg<<endl;
    return count;
}
