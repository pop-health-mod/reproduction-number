# reproduction-number

Estimation du taux de reproduction du SRAS-CoV-2 au Québec —

# Précisions méthodologiques
## Qu’est que le taux de reproduction effectif (R\_t)?

Le taux de reproduction effectif (R\_t) est défini comme le nombre moyen d’infections secondaires produites par un cas dans une population où certains individus ne sont plus susceptibles à l’infection. 

## Comment interpréter R\_t et pourquoi le calculer ?

Un R\_t plus grand que 1 indique une accélération exponentielle de la transmission, tandis qu’un R\_t plus petit que 1 indique plutôt que l’épidémie est en recul. Lorsque le R\_t est près de la valeur unitaire, la transmission et le nombre de cas demeurent stable. Le R\_t est une mesure qui permet de détecter les changements du niveau de transmission d’un agent infectieux dans le temps. L’objectif est de maintenir ce taux de reproduction en dessous de 1 en adaptant et optimisant les interventions pour réduire la propagation du SRAS-CoV-2.
Il faut interpréter cet indicateur en fonction du contexte de transmission. À titre d’exemple, si l’infection est bien contrôlée (c.-à-d. peu de cas) et qu’on observe une éclosion, le R\_t augmentera. Une valeur élevée de R\_t est cependant beaucoup plus inquiétante si le nombre de cas est déjà lui-même déjà élevé. Notre interprétation du R\_t doit donc tenir compte des autres données épidémiologiques disponibles et son interprétation nous donne une information complémentaire.

## Sommaire méthodologique

La méthodologie adoptée permet d’estimer un R\_t dit « instantané » (1,2), à partir des cas de SRAS-CoV-2 détectés au Québec. Cette méthodologie requiert deux types de données :

 - La courbe épidémiologique des nouvelles infections.
 - L’intervalle générationnel (le délai moyen entre l’infection d’un cas primaire et celle d’un cas secondaire; Figure 1).

Comme il est impossible d’observer directement le premier type de données – la courbe des infections (c.-à-d., la date à laquelle une personne est infectée) – elle est calculée rétroactivement à partir de la série temporelle des cas confirmés en laboratoire (date du prélèvement) ou par lien épidémiologique (date de déclaration). Cette série temporelle est d’abord lissée afin de réduire l’impact des effets de fins de semaine où un nombre de cas plus faible est généralement déclaré (3). La série lissée est ensuite utilisée pour calculer la courbe des infections à partir d’un algorithme de déconvolution de type Richardson-Lucy (4) qui permet d’estimer les temps d’infection. Cet algorithme s’appuie sur le délai entre l’infection et l’apparition des symptômes (période d’incubation; Figure 1) ainsi que les délais entre l’apparition des symptômes et le prélèvement du test (ou la date de déclaration des cas par lien épidémiologique). Cette méthode permet également de tenir compte de la troncature à droite des séries temporelles (c.-à-d., les personnes infectées, mais pas encore confirmées) (4,5). Le R\_t est ensuite calculé, sur une fenêtre mobile de 5 jours, comme le ratio des nouvelles infections sur le profil d’infectiosité dans la population au même moment (1,2). Le profil d’infectiosité est une distribution qui représente la probabilité journalière qu’un individu soit infectieux avoir été infecté. À noter que les cas importés contribuent au profil d’infectiosité mais sont exclus du numérateur de ce ratio (ne sont pas considérés comme de « nouvelles » infections). Il existe des méthodes alternatives pour calculer le taux de reproduction mais la méthode du R\_t « instantané » est la plus robuste pour obtenir des estimés en temps quasi-réel (5).

## Limites méthodologiques

Estimer le R\_t  de façon précise présente plusieurs défis (5). Premièrement, les efforts de dépistage et les critères d’éligibilité au test du SRAS-CoV-2 peuvent varier dans le temps, influençant le nombre de cas diagnostiqués. Une augmentation ou une diminution importante des efforts de dépistage pourrait causer un accroissement ou une réduction artificielle du R\_t, respectivement, puisque le modèle interprète tous changements dans le nombre de cas comme une modification du taux de transmission. Deuxièmement, les délais d’agrégation des données de vigie épidémiologique nous obligent à exclure les deux dernières journées des séries temporelles. Ceux-ci se rajoutent aux délais épidémiologiques (incubation et prélèvement) faisant en sorte que le R\_t estimé représente le niveau de transmission observé il y a près de 10 jours. Troisièmement, certaines des distributions utilisées (intervalle générationnel, symptômes-prélèvement) pour reconstruire les courbes des infections pourraient varier dans le temps et affecter les estimés de R\_t. Si nécessaire, des itérations futures nous permettront de prendre en compte ces changements potentiels. Finalement, lorsque le niveau de transmission est bas, le 𝓡t est moins utile comme indicateur épidémiologique. Le cas échéant, la proportion des cas pouvant être liée à une éclosion connue est un meilleur indicateur du niveau de transmission.

## Télécharger le repo




## Équipe de modélisation
Mathieu Maheu-Giroux
Arnaud Godin
Yiqing Xia

## Collaborateurs
Marc Brisson (Université Laval)
David Buckeridge (Université McGill)
Charlotte Lanièce (Université McGill)
Marie-Claude Boily (Imperial College London)

