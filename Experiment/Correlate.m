%% INPUT: a sequence of rx data
%% Type: 1 short, 2 long, 3 cp+long, 4 rawofdm, 0 self-defined seq
function conv_res = Correlate(rxdata, type, seq)
  
  % 11a style
  short = [
   0.0460 + 0.0460i;
  -0.1324 + 0.0023i;
  -0.0135 - 0.0785i;
   0.1428 - 0.0127i;
   0.0920          ;
   0.1428 - 0.0127i;
  -0.0135 - 0.0785i;
  -0.1324 + 0.0023i;
   0.0460 + 0.0460i;
   0.0023 - 0.1324i;
  -0.0785 - 0.0135i;
  -0.0127 + 0.1428i;
        0 + 0.0920i;
  -0.0127 + 0.1428i;
  -0.0785 - 0.0135i;
   0.0023 - 0.1324i
   ];

  cp = [
   -0.1562          ;
   0.0123 - 0.0976i;
   0.0917 - 0.1059i;
  -0.0919 - 0.1151i;
  -0.0028 - 0.0538i;
   0.0751 + 0.0740i;
  -0.1273 + 0.0205i;
  -0.1219 + 0.0166i;
  -0.0350 + 0.1509i;
  -0.0565 + 0.0218i;
  -0.0603 - 0.0813i;
   0.0696 - 0.0141i;
   0.0822 - 0.0924i;
  -0.1313 - 0.0652i;
  -0.0572 - 0.0393i;
   0.0369 - 0.0983i;
   0.0625 + 0.0625i;
   0.1192 + 0.0041i;
  -0.0225 - 0.1607i;
   0.0587 + 0.0149i;
   0.0245 + 0.0585i;
  -0.1368 + 0.0474i;
   0.0010 + 0.1150i;
   0.0533 - 0.0041i;
   0.0975 + 0.0259i;
  -0.0383 + 0.1062i;
  -0.1151 + 0.0552i;
   0.0598 + 0.0877i;
   0.0211 - 0.0279i;
   0.0968 - 0.0828i;
   0.0397 + 0.1112i;
  -0.0051 + 0.1203i
  ];

  long = [
    0.1562         ;
  -0.0051 - 0.1203i;
   0.0397 - 0.1112i;
   0.0968 + 0.0828i;
   0.0211 + 0.0279i;
   0.0598 - 0.0877i;
  -0.1151 - 0.0552i;
  -0.0383 - 0.1062i;
   0.0975 - 0.0259i;
   0.0533 + 0.0041i;
   0.0010 - 0.1150i;
  -0.1368 - 0.0474i;
   0.0245 - 0.0585i;
   0.0587 - 0.0149i;
  -0.0225 + 0.1607i;
   0.1192 - 0.0041i;
   0.0625 - 0.0625i;
   0.0369 + 0.0983i;
  -0.0572 + 0.0393i;
  -0.1313 + 0.0652i;
   0.0822 + 0.0924i;
   0.0696 + 0.0141i;
  -0.0603 + 0.0813i;
  -0.0565 - 0.0218i;
  -0.0350 - 0.1509i;
  -0.1219 - 0.0166i;
  -0.1273 - 0.0205i;
   0.0751 - 0.0740i;
  -0.0028 + 0.0538i;
  -0.0919 + 0.1151i;
   0.0917 + 0.1059i;
   0.0123 + 0.0976i;
  -0.1562          ;
   0.0123 - 0.0976i;
   0.0917 - 0.1059i;
  -0.0919 - 0.1151i;
  -0.0028 - 0.0538i;
   0.0751 + 0.0740i;
  -0.1273 + 0.0205i;
  -0.1219 + 0.0166i;
  -0.0350 + 0.1509i;
  -0.0565 + 0.0218i;
  -0.0603 - 0.0813i;
   0.0696 - 0.0141i;
   0.0822 - 0.0924i;
  -0.1313 - 0.0652i;
  -0.0572 - 0.0393i;
   0.0369 - 0.0983i;
   0.0625 + 0.0625i;
   0.1192 + 0.0041i;
  -0.0225 - 0.1607i;
   0.0587 + 0.0149i;
   0.0245 + 0.0585i;
  -0.1368 + 0.0474i;
   0.0010 + 0.1150i;
   0.0533 - 0.0041i;
   0.0975 + 0.0259i;
  -0.0383 + 0.1062i;
  -0.1151 + 0.0552i;
   0.0598 + 0.0877i;
   0.0211 - 0.0279i;
   0.0968 - 0.0828i;
   0.0397 + 0.1112i;
  -0.0051 + 0.1203i
  ];

  rawofdm_sts_a = [
   1.2984 - 0.8548j;
  -1.0196 + 0.3308j;
  -0.1185 + 0.3173j; 
   0.5303 + 0.8118j;
  -0.1019 + 0.7091j;
   0.2147 - 0.4824j;
  -0.6892 + 0.3757j;
   0.4383 + 0.0783j;
   0.4786 - 1.1558j;
  -0.1643 - 0.4001j;
   0.3664 + 0.6463j;
  -0.3522 + 0.4780j;
   0.2358 + 0.0465j;
   0.4488 + 0.5154j;
   0.0414 + 1.0010j;
   0.1826 + 1.5399j;
  -0.1926 + 0.9062j;
  -0.1387 - 0.8931j;
  -0.1599 - 0.6938j;
   0.1183 - 0.2761j;
   0.0605 - 0.3776j;
   0.3285 + 0.2880j;
   0.5988 - 0.1381j;
  -0.5900 - 1.3076j;
  -0.7603 - 1.0833j;
  -0.7019 + 0.0782j;
  -0.4639 - 0.0143j;
  -0.3260 - 0.3310j;
   0.5487 - 0.2704j;
   1.6965 - 0.3181j;
  -1.1422 + 0.5860j;
  -0.6653 - 0.1119j
  ];

  rawofdm_sts_b = [
 -2.670712947845459 - 2.180722713470459i
 -0.651309609413147 + 1.010068058967590i
  0.527148604393005 + 1.763299226760864i
  1.970795154571533 - 0.972449064254761i
 -1.084294557571411 - 0.210654646158218i
 -1.231719493865967 + 0.594141066074371i
 -0.646610617637634 - 1.151523709297180i
 -1.194927215576172 - 0.873271703720093i
  0.479463338851929 - 1.442418932914734i
  0.283250272274017 - 0.732088685035706i
 -0.014967918395996 + 0.870575070381165i
 -0.883156239986420 - 1.091405987739563i
  0.841613650321960 + 1.558691740036011i
  1.162610769271851 + 2.308510541915894i
 -2.203146219253540 - 1.207847118377686i
  0.550024509429932 + 1.909584999084473i
  0.973607599735260 + 2.049424171447754i
 -1.049502491950989 - 0.651550054550171i
  1.411922574043274 - 0.662048041820526i
  2.327744960784912 - 2.089462041854858i
  0.974457740783691 - 0.650885224342346i
  0.304966032505035 + 0.454963862895966i
  2.090454101562500 - 0.684922099113464i
  0.659242093563080 - 0.879358410835266i
 -0.866571307182312 - 0.645443797111511i
  0.273810207843781 + 1.644728302955627i
 -1.712809562683105 + 1.624717116355896i
  0.456649601459503 - 0.386361002922058i
  1.674814581871033 - 0.070392310619354i
 -0.167259514331818 + 0.682005345821381i
  0.225631237030029 + 1.040151357650757i
 -2.811218738555908 - 0.928055703639984i
  ];
  rawofdm_sts_b = rawofdm_sts_b ./ 2;

  rawofdm_lts_a = [
 -0.324729651212692 - 0.648506164550781i
 -0.301016688346863 + 0.359643161296844i
 -0.147326365113258 + 1.058500766754150i
  0.308397918939590 + 0.353920727968216i
  0.113463878631592 - 0.555182576179504i
  0.783052086830139 - 0.210373610258102i
  1.038782119750977 + 0.523121595382690i
  0.038778871297836 + 1.057318687438965i
  0.303817242383957 + 0.611358046531677i
  0.351564615964890 + 0.193092763423920i
  0.042094588279724 + 0.011735804378986i
 -0.725173592567444 - 0.853282272815704i
 -1.025235891342163 + 0.565869688987732i
 -0.180808752775192 + 1.309297442436218i
 -0.232895791530609 - 0.031507849693298i
  0.803657352924347 + 0.023751318454742i
  0.809926629066467 + 0.388462454080582i
 -0.696269154548645 + 1.213965773582458i
  0.536057114601135 - 0.360280871391296i
  0.534684121608734 - 1.910227775573730i
 -0.368418514728546 + 0.425022453069687i
  0.956639289855957 + 0.734396398067474i
  0.715041041374207 + 0.537298798561096i
 -0.114415407180786 + 0.320772647857666i
  0.258616387844086 - 0.332280337810516i
  0.405675143003464 + 0.680152058601379i
  0.968100428581238 - 0.588708758354187i
  0.505425691604614 - 0.558997452259064i
 -0.708929419517517 + 0.109176024794579i
 -0.107907086610794 - 1.016447663307190i
  0.014984667301178 + 0.613836824893951i
 -0.272099196910858 + 0.367176473140717i
  0.267610937356949 - 1.178259372711182i
 -0.502658367156982 - 0.648753941059113i
 -0.524189591407776 - 0.665371358394623i
  0.179502278566360 + 0.150814265012741i
 -0.793793082237244 + 0.505022168159485i
 -1.009974360466003 - 0.379487127065659i
  0.107874453067780 - 0.148828387260437i
  0.290660738945007 + 0.710338413715363i
 -0.821817219257355 + 0.384588837623596i
  0.023147732019424 - 0.065792113542557i
  1.620932698249817 + 0.389522135257721i
 -0.228508919477463 - 0.353674590587616i
 -0.821103811264038 - 0.008949890732765i
  0.688975334167480 + 0.764788985252380i
  0.194231986999512 + 0.098834812641144i
 -0.246814742684364 + 0.934092581272125i
 -0.280057817697525 + 0.248521834611893i
 -1.127288818359375 - 1.080757260322571i
 -1.408289909362793 - 0.583601474761963i
 -0.066827774047852 - 0.672649919986725i
 -0.003877311944962 - 0.231080800294876i
 -0.951217770576477 - 0.248973965644836i
  0.549793362617493 - 0.196395725011826i
  0.285471946001053 + 0.495385587215424i
 -0.476266354322433 - 0.027247756719589i
  0.536776483058929 - 0.045307524502277i
 -0.651929199695587 + 0.185041368007660i
 -0.215558499097824 + 0.494138538837433i
  0.370194315910339 + 0.114210337400436i
 -0.757714629173279 - 1.359393239021301i
  0.677338123321533 - 1.213922619819641i
  0.811844050884247 - 0.758927166461945i];

  rawofdm_lts_b = [
  0.600603103637695 - 0.642856299877167i
  0.555812895298004 + 0.004549294710159i
  0.584578216075897 + 0.615208506584167i
 -0.158044815063477 + 0.712308406829834i
 -0.960626363754272 - 0.045826829969883i
 -0.660843372344971 + 0.553179442882538i
 -0.514051437377930 + 0.805295705795288i
 -0.470963001251221 + 0.542817950248718i
 -0.157838299870491 + 0.367176592350006i
  0.361174076795578 - 0.239544808864594i
  0.132994443178177 + 0.043339848518372i
 -0.630867958068848 + 0.092171967029572i
 -0.112968862056732 + 0.529931843280792i
 -0.312588065862656 + 0.876693308353424i
 -0.275063216686249 + 0.800422430038452i
  0.361452966928482 + 0.974227130413055i
 -0.522131204605103 - 0.099471867084503i
 -0.126634806394577 + 0.235434323549271i
  0.062711864709854 + 0.397178828716278i
 -0.441921979188919 - 0.809287309646606i
  0.653882503509521 - 0.261922299861908i
  0.915526270866394 - 0.259624749422073i
  0.373903781175613 - 0.123748466372490i
 -0.182362228631973 - 0.272522628307343i
 -0.247716516256332 - 1.031872987747192i
  0.130689084529877 + 0.174283951520920i
 -0.088178709149361 - 0.294190526008606i
 -0.395977854728699 - 0.813577950000763i
 -0.007834650576115 + 0.957274079322815i
  1.360511064529419 + 0.837775111198425i
  0.145312249660492 - 0.354911625385284i
 -0.910548269748688 - 0.213123261928558i
  0.456434339284897 + 0.589231312274933i
 -1.513518571853638 + 0.267002433538437i
 -1.148481369018555 + 0.013660699129105i
  1.170813083648682 + 1.133775949478149i
 -0.464080572128296 + 0.128193050622940i
  0.086121529340744 - 1.316248655319214i
  0.523923635482788 + 0.232003122568130i
 -0.518945693969727 + 0.118401855230331i
  0.494644552469254 - 1.870595455169678i
  0.539611577987671 - 0.933428049087524i
  0.117953240871429 + 0.682337582111359i
  0.012290954589844 + 1.330901145935059i
  0.113877646625042 + 0.735423684120178i
  0.237517833709717 - 1.392126679420471i
  0.063391976058483 - 0.318477690219879i
  0.254892766475677 + 1.169041395187378i
 -0.226137489080429 - 0.086121782660484i
  0.204377382993698 - 0.703465938568115i
  0.658793449401855 - 1.584909915924072i
  0.095155924558640 - 0.721519708633423i
  0.243911772966385 + 0.322156012058258i
 -0.281243562698364 - 0.556648969650269i
 -0.430657893419266 + 0.189154192805290i
 -0.392031759023666 - 0.578258812427521i
 -0.508833467960358 - 1.006151914596558i
  0.278049200773239 - 0.324500679969788i
 -0.308683604001999 - 0.994299888610840i
 -0.341718554496765 + 1.144831657409668i
  0.683588445186615 + 1.349358081817627i
  0.249712139368057 - 0.672454535961151i
  0.062778450548649 - 0.041987776756287i
  0.524501562118530 - 0.361062824726105i];

  if type == 1
      signal = short;
  elseif type == 2
      signal = long;
  elseif type == 3
      signal = [cp; long];
  elseif type == 4
      signal = rawofdm_sts_a;
  elseif type == 5
      signal = rawofdm_sts_b;
  elseif type == 6
      signal = rawofdm_lts_a;
  elseif type == 7
      signal = rawofdm_lts_b;
  elseif type == 0
      signal = seq;
  end
  
  conv_res = filter(conj(signal(end:-1:1)), 1, rxdata);
  conv_res = abs(conv_res);
%   figure;
%   plot(conv_res,'b.-');

end