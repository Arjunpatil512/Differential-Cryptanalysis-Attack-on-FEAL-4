/** ======================================================================================
    EXPERIMENT: Criptoanàlisi diferencial del criptosistema FEAL-4.
    AUTOR:      Eric Santiño Cervera

    EXPLICACIÓ: Utilitza l'atac diferencial explicat durant el treball per a
                aconseguir recuperar les 6 subclaus emprades en el procés de
                xifratge (KL, KR, K1, K2, K3, K4).
======================================================================================= */
//-----------------------------------------------------------------
// NOTACIONS:   P  = plaintext / missatge pla
//              C  = ciphertext / missatge xifrat
//              Kj = claus de ronda j
//
//              L  = lsub-bloc esquerre
//              R  = sub-bloc dret
//
//              dP = diferència de missatges plans
//              dC = diferència de missatges xifrats
//              dX = diferència d'entrada
//              dY = diferència de sortida
//-----------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NUM_PARELLES 8

typedef unsigned char uint8;
typedef unsigned long uint32;
typedef unsigned long long uint64;

//==================================
//   VARIABLES GLOBALS
//==================================
// Claus de ronda.
uint32 K[4], KL, KR;

// Missatges plans i xifrats per a trencar el criptosistema.
uint64 P_0[NUM_PARELLES], P_1[NUM_PARELLES];
uint64 C_0[NUM_PARELLES], C_1[NUM_PARELLES];

//==================================
//   FUNCIONS AUXILIARS
//==================================
// Manipulació dels blocs de 64 bits, i els sub-blocs que els composen.
uint32 diff_MeitatEsquerra(uint64 bloc) {
    return (uint32)(bloc >> 32);
}
uint32 diff_MeitatDreta(uint64 bloc) {
    return (uint32)(bloc & 0xFFFFFFFF);
}
uint64 diff_CombinarMeitats(uint32 meitat_esquerra, uint32 meitat_dreta) {
    return ((uint64)(meitat_esquerra) << 32) | meitat_dreta;
}

// Manipulació dels sub-blocs de 32 bits, i els bytes que els composen.
uint8 diff_SepararBytes(uint32 x, uint8 n) {
    return (uint8)((x >> (8*n)) & 0xFF);
}
uint32 diff_CombinarBytes(uint8 *bytes) {
    uint32 temp, temp2;

    temp  = ((uint32)(bytes[3]) << 24) | ((uint32)(bytes[2]) << 16);
    temp2 = ((uint32)(bytes[1]) << 8 ) |  (uint32)(bytes[0]);

    return temp | temp2;
}

// Funció de ronda de la xarxa de Feistel.
uint8 diff_func_G(uint8 a, uint8 b, uint8 x) {
    x = (uint8)((a + b + x) % 256);
    return (x << 2) | (x >> 6);
}
uint32 diff_func_f(uint32 R) {
    uint8 k, x[4], y[4];

    for(k=0; k<4; k++) {
        x[k] = diff_SepararBytes(R, k);
    }

    y[1] = diff_func_G(x[0]^x[1], x[2]^x[3], 1);
    y[0] = diff_func_G(x[0]     , y[1]     , 0);
    y[2] = diff_func_G(x[2]^x[3], y[1]     , 0);
    y[3] = diff_func_G(x[3]     , y[2]     , 1);

    return diff_CombinarBytes(y);
}

// Procesos de xifratge i desxifratge del FEAL-4.
void diff_XOR_dreta(uint32 *L, uint32 *R) {
    *R = *L ^ *R;
}
uint64 diff_Xifrar(uint64 P) {
    uint8 i;
    uint32 L, R, nou_L, nou_R;

    L = diff_MeitatEsquerra(P) ^ KL;
    R = diff_MeitatDreta(P) ^ KR;
    diff_XOR_dreta(&L, &R);

    for(i=0; i<4; i++) {
        if (i!=3) {
            nou_L = R;
            nou_R = L ^ diff_func_f(R ^ K[i]);
        }
        else {
            nou_L = L ^ diff_func_f(R ^ K[i]);
            nou_R = R;
        }

        L = nou_L; R = nou_R;
    }

    diff_XOR_dreta(&L, &R);
    return diff_CombinarMeitats(L, R);
}
void diff_DesxifrarRonda(uint64 *missatge, uint8 index) {
    uint32 L, R;

    L = diff_MeitatDreta(*missatge);
    R = diff_MeitatEsquerra(*missatge) ^ diff_func_f(L ^ K[index-1]);

    *missatge = diff_CombinarMeitats(L, R);
}

// Generació de valors aleatoris: claus de ronda, parelles de missatges...
void diff_GenerarClaus(uint32 llavor) {
    uint8 i;
    srand(llavor);

    for(i=0; i<4; i++) {
        K[i] = ((uint32)(rand()) << 16) | ((uint32)(rand()) & 0xFFFF);
    }

    KL = ((uint32)(rand()) << 16) | ((uint32)(rand()) & 0xFFFF);
    KR = ((uint32)(rand()) << 16) | ((uint32)(rand()) & 0xFFFF);
}
void diff_GenerarParelles(uint64 dP) {
    uint8 k;

    for(k=0; k<NUM_PARELLES; k++) {
        P_0[k]  = (uint64)(rand() & 0xFFFF) << 48;
        P_0[k] ^= (uint64)(rand() & 0xFFFF) << 32;
        P_0[k] ^= (uint64)(rand() & 0xFFFF) << 16;
        P_0[k] ^= (uint64)(rand() & 0xFFFF);
        P_1[k]  = P_0[k] ^ dP;

        C_0[k] = diff_Xifrar(P_0[k]);
        C_1[k] = diff_Xifrar(P_1[k]);
    }
}

// Atac diferencial per obtenir la última clau utilitzada.
uint64 diff_AtacDiferencial(uint32 dXrondaL) {
    uint8 k;
    uint32 candidat_H, contador_H;
    uint32 L_0, L_1, R_0, R_1, fronda_0, fronda_1;

    for(candidat_H=0x00000000L; candidat_H<0xFFFFFFFFL; candidat_H++) {
        contador_H=0;

        for(k=0; k<NUM_PARELLES; k++) {
            L_0 = diff_MeitatEsquerra(C_0[k]); R_0 = diff_MeitatDreta(C_0[k]);
            L_1 = diff_MeitatEsquerra(C_1[k]); R_1 = diff_MeitatDreta(C_1[k]);

            fronda_0 = diff_func_f(R_0 ^ candidat_H);
            fronda_1 = diff_func_f(R_1 ^ candidat_H);

            if (dXrondaL == (L_0 ^ fronda_0 ^ L_1 ^ fronda_1)) {
                contador_H++;
            }
            else break;
        }
        if (contador_H == NUM_PARELLES) {
            break;
        }
    }

    return candidat_H;
}

// =========================================
// ATAC DIFERENCIAL COMPLERT DEL FEAL-4.
// =========================================
int main() {
    uint8 i, j, k;
    uint32 clau_trobada, candidat_K1, candidat_KL, candidat_KR;
    uint32 temps_inicial, temps_inicial_ronda, L_0, L_1, R_0, R_1, dX;
    uint64 dP[] = {0x8080000080800000LL, 0x0000000080800000LL, 0x0000000002000000LL};

    printf("================================================================\n");
    printf("   *** CRIPTOANALISI DIFERENCIAL DEL CRIPTOSISTEMA FEAL-4 ***  \n");
    printf("================================================================\n");
    diff_GenerarClaus(time(NULL));
    temps_inicial = time(NULL);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ATAC A LES RONDES 4, 3, 2
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for(i=4; i>=2; i--) {
        printf("-> RONDA %i\n", i);
        temps_inicial_ronda = time(NULL);

        // Generem les parelles de missatges plans i xifrats.
        diff_GenerarParelles(dP[4-i]);
        printf("   Generant %i parelles de missatges amb:\n", NUM_PARELLES);
        printf("   dX1 = 0x%016llX\n   dY%i = 0x02000000********\n\n", dP[4-i], i-1);

        // Desfem el canvi final, i tantes rondes com claus de ronda ja hem extret.
        for(k=0; k<NUM_PARELLES; k++) {
            L_0 = diff_MeitatEsquerra(C_0[k]); R_0 = diff_MeitatDreta(C_0[k]);
            L_1 = diff_MeitatEsquerra(C_1[k]); R_1 = diff_MeitatDreta(C_1[k]);
            diff_XOR_dreta(&L_0, &R_0);
            diff_XOR_dreta(&L_1, &R_1);

            C_0[k] = diff_CombinarMeitats(L_0, R_0);
            C_1[k] = diff_CombinarMeitats(L_1, R_1);

            for(j=4; j>i; j--) {
                diff_DesxifrarRonda(&C_0[k], j);
                diff_DesxifrarRonda(&C_1[k], j);
            }
        }

        // Trobem la clau d'aquesta ronda.
        clau_trobada = diff_AtacDiferencial(0x02000000L);

        if (clau_trobada == K[i-1]) {
            printf("     Clau: 0x%08lX", clau_trobada);
            printf("   | EXIT!   Temps d'execucio per a trobar la clau de ronda K%i: %i segons\n\n", i, time(NULL) - temps_inicial_ronda);
        }
        else {
            printf("     Clau: 0x%08lX", clau_trobada);
            printf("   | FRACAS! No s'ha pogut trobar la clau de ronda K%i", i);
            return -1;
        }
    }

    //~~~~~~~~~~~~~~~~~~~~~~
    // ATAC A LA RONDA 1
    //~~~~~~~~~~~~~~~~~~~~~~
    printf("-> RONDA 1\n");
 	printf("   Provant valors de claus de ronda fins a obtenir un que funcioni\n");
	printf("   per a tots els parells de missatges que tenim...\n\n");
    temps_inicial_ronda = time(NULL);

 	for(k=0; k<NUM_PARELLES; k++) {
        diff_DesxifrarRonda(&C_0[k], 2);
        diff_DesxifrarRonda(&C_1[k], 2);
 	}

    for(candidat_K1=0; candidat_K1<0xFFFFFFFFL; candidat_K1++) {
        candidat_KL = 0;
        candidat_KR = 0;

        // Volem trobar ara el valor de les claus que ens passin de C_0 a P_0,
        // per a tots els parells de missatges que tenim.
        for(k=0; k<NUM_PARELLES; k++) {

            // 0 és ara el missatge pla P_0.
            L_0 = diff_MeitatEsquerra(P_0[k]);
            R_0 = diff_MeitatDreta(P_0[k]);

            // 1 és ara el missatge xifrat C_0.
            L_1 = diff_MeitatEsquerra(C_0[k]);
            R_1 = diff_MeitatDreta(C_0[k]);

	 	   	dX = diff_func_f(R_1 ^ candidat_K1) ^ R_0;

            if (candidat_KL == 0) {
                candidat_KL = dX ^ L_0;
                candidat_KR = dX ^ R_1 ^ R_0;
            }
            else if (((dX ^ L_0) != candidat_KL) || ((dX ^ R_1 ^ R_0) != candidat_KR)) {
                break;
            }
        }
        if ((candidat_KL == KL) && (candidat_KR == KR)) {
            printf("     Clau K1: 0x%08lX\n", candidat_K1);
            printf("     Clau KL: 0x%08lX\n", candidat_KL);
            printf("     Clau KR: 0x%08lX", candidat_KR);
            printf(" | EXIT!   Temps d'execucio per a trobar les claus de ronda: %i segons\n\n", time(NULL) - temps_inicial_ronda);
            break;
        }
    }

    printf("** Criptoanalisi diferencial finalitzada amb exit! **\n");
    printf("      KL = 0x%08LX\n", KL);
    printf("      KR = 0x%08LX\n", KR);
    printf("      K1 = 0x%08LX\n", K[0]);
    printf("      K2 = 0x%08LX\n", K[1]);
    printf("      K3 = 0x%08LX\n", K[2]);
    printf("      K4 = 0x%08LX\n\n", K[3]);
    printf("   TEMPS TOTAL D'EXECUCIO = %i segons\n\n", time(NULL) - temps_inicial);
    return 0;
}

/** -----------------------------------------------------------------------------------------
        Idea original per Jon King ( http://theamazingking.com )
        Adaptat i re-escrit
------------------------------------------------------------------------------------------ */
