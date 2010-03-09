%%
%% This is file `achemso.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% achemso.dtx  (with options: `class')
%%  ----------------------------------------------------------------
%%  achemso --- Support for submissions to American  Chemical
%%    Society journals
%%  Maintained by Joseph Wright
%%  E-mail: joseph.wright@morningstar2.co.uk
%%  Released under the LaTeX Project Public License v1.3c or later
%%  See http://www.latex-project.org/lppl.txt
%%  ----------------------------------------------------------------
%% 
\NeedsTeXFormat{LaTeX2e}
\LoadClass[12pt]{article}
\RequirePackage[etex=false]{notes2bib}[2008/06/21]
\RequirePackage{achemso}
\ProvidesClass{achemso}
  [\acs@ver Submissions to ACS journals]
\RequirePackage[T1]{fontenc}
\RequirePackage[scaled=0.90]{helvet}
\RequirePackage[margin=2.54cm]{geometry}
\RequirePackage{mathptmx,courier,setspace,graphicx,truncate,%
  float,varioref}
\AtBeginDocument{\doublespacing}
\renewcommand*{\title}[2][]{%
  \gdef\@title{#2}%
  \gdef\acs@title{#2}%
  \gdef\@shorttitle{#1}}
\@onlypreamble\title
\newcount\acs@authorcnt
\newcount\acs@affilcnt
\newcount\acs@altaffilcnt
\acs@affilcnt\@ne\relax
\acs@altaffilcnt\@ne\relax
\renewcommand*{\author}[1]{%
  \global\advance\acs@authorcnt\@ne\relax
  \expandafter\gdef
    \csname @author@\@roman\the\acs@authorcnt\endcsname{#1}%
  \begingroup
    \advance\acs@affilcnt\@ne\relax
    \expandafter\xdef
      \csname @author@affil@\@roman\the\acs@authorcnt\endcsname
        {\the\acs@affilcnt}%
  \endgroup}
\@onlypreamble\author
\renewcommand*{\and}{%
  \ClassError{achemso}{\string\and\space not supported}
    {The achemso class does not use \string\and\MessageBreak
     see the documentation for details}}
\renewcommand*{\thanks}[1]{%
  \ClassError{achemso}{\string\thanks\space not supported}
    {The achemso class does not use \string\thanks\MessageBreak
     see the documentation for details}}
\newcommand*{\affiliation}[2][\relax]{%
  \ifnum\acs@authorcnt>\z@\relax
    \global\advance\acs@affilcnt\@ne
    \begingroup
      \expandafter\def
        \csname @address@\@roman\acs@affilcnt\endcsname{#2}%
      \acs@tempcnta\acs@affilcnt\relax
      \acs@ifdupaffil
        {\begingroup
           \acs@tempcntb\@ne\relax
           \acs@switchfalse
           \edef\acs@tempa{%
             \csname @address@\@roman\acs@tempcnta\endcsname}%
           \acs@ifdup@affil
           \acs@tempcnta\z@\relax
           \edef\acs@tempa{\the\acs@affilcnt}%
           \global\advance\acs@affilcnt\m@ne\relax
           \acs@swapaffil
         \endgroup}
        {\expandafter\gdef
           \csname @address@\@roman\acs@affilcnt\endcsname{#2}%
         \ifx\relax#1\relax
           \expandafter\gdef
             \csname @affil@\@roman\acs@affilcnt\endcsname{#2}%
         \else
           \expandafter\gdef
             \csname @affil@\@roman\acs@affilcnt\endcsname{#1}%
         \fi}
    \endgroup
  \else
    \ClassWarning{achemso}
      {Affiliation with no author}%
  \fi}
\@onlypreamble\affiliation
\newcommand*{\acs@swapaffil}{%
  \advance\acs@tempcnta\@ne\relax
  \ifnum\acs@tempcnta>\acs@authorcnt\relax\else
    \edef\acs@tempb{%
      \csname @author@affil@\@roman\acs@tempcnta\endcsname}%
    \ifx\acs@tempa\acs@tempb
      \expandafter\xdef
        \csname @author@affil@\@roman\acs@tempcnta\endcsname{%
          \the\acs@tempcntb}%
    \fi
    \acs@swapaffil
  \fi}
\newcommand*{\altaffiliation}[1]{%
  \ifnum\acs@authorcnt>\z@\relax
    \global\advance\acs@altaffilcnt\@ne\relax
    \expandafter\gdef
      \csname @altaffil@\@roman\acs@authorcnt\endcsname{#1}%
    \expandafter\xdef
      \csname @author@altaffil@\@roman\acs@authorcnt\endcsname
        {\the\acs@altaffilcnt}
  \else
    \ClassWarning{achemso}
      {Affiliation with no author}%
  \fi}
\@onlypreamble\altaffiliation
\newcommand*{\email}[1]{%
  \ifnum\acs@authorcnt>\z@\relax
    \expandafter\gdef
      \csname @email@\@roman\acs@authorcnt\endcsname{#1}%
  \else
    \ClassWarning{achemso}
      {E-mail with no author}%
  \fi}
\@onlypreamble\email
\renewcommand*{\@maketitle}{%
  \ifnum\acs@authorcnt<\z@\relax
    \ClassError{achemso}{No authors defined}
      {At least one author is required}%
  \else
    \newpage
    \null
    \vskip 2em%
    \begin{center}%
      {\LARGE\bfseries\sffamily
       \renewcommand*{\acs@tempa}{suppinfo}%
       \ifx\acs@manuscript\acs@tempa
         Supporting information for:
       \fi
       \@title \par}%
      \vskip 1.5em\relax
      {\large\sffamily\frenchspacing \acs@authorlist}%
      \vskip 1em%
      {\itshape\acs@addresslist}%
      \ifnum\acs@affilcnt>\tw@\relax
        \acs@affilfoot
      \else
        \ifnum\acs@altaffilcnt>\@ne\relax
          \acs@affilfoot
        \fi
      \fi
      \vskip 1em\relax
      {\sffamily E-mail: \acs@emaillist}%
    \end{center}
    \par
    \vskip 1.5em\relax
  \fi}
\newcommand*{\acs@authorlist}{%
  \acs@tempcnta\z@\relax
  \acs@tempcntc\z@\relax
  \acs@author@list}
\newcommand*{\acs@author@list}{%
  \advance\acs@tempcnta\@ne\relax
  \ifnum\acs@tempcnta>\acs@authorcnt\relax\else
    \ifnum\acs@tempcnta=\acs@authorcnt\relax
      \ifnum\acs@tempcnta=\@ne\relax\else
        and
      \fi
    \fi
    \csname @author@\@roman\acs@tempcnta\endcsname
    \ifnum\acs@tempcnta=\acs@authorcnt\relax\else
      ,%
    \fi
    \begingroup
      \@ifundefined{@email@\@roman\acs@tempcnta}
        {\aftergroup\@firstoftwo}
        {\aftergroup\@secondoftwo}%
    \endgroup
      {\def\acs@tempb{}}
      {\protected@edef\acs@tempb{%
         \acs@fnsymbol{\@ne}%
         \ifnum\acs@affilcnt>\tw@\relax
           ,%
         \else
           \ifnum\acs@altaffilcnt>\@ne\relax
           ,%
           \fi
         \fi}}%
    \ifnum\acs@affilcnt>\tw@\relax
      \protected@edef\acs@tempb{\acs@tempb\@fnsymbol{%
        \csname @author@affil@\@roman\acs@tempcnta
          \endcsname}}%
    \else
      \ifnum\acs@altaffilcnt>\@ne\relax
        \protected@edef\acs@tempb{\acs@tempb\@fnsymbol{%
          \csname @author@affil@\@roman\acs@tempcnta
            \endcsname}}%
      \fi
    \fi
    \begingroup
      \@ifundefined{@author@altaffil@\@roman\acs@tempcnta}
        {\aftergroup\@gobble}
        {\aftergroup\@firstofone}%
    \endgroup
      {\global\advance\acs@tempcntc\@ne\relax
       \advance\acs@tempcntc\acs@affilcnt
       \ifnum\acs@affilcnt>\@ne\relax
         \protected@edef\acs@tempb{\acs@tempb,}%
       \fi
       \protected@edef\acs@tempb{%
         \acs@tempb\@fnsymbol{\acs@tempcntc}}}%
    \textsuperscript{\acs@tempb}
    \acs@author@list
  \fi}
\newcommand*{\acs@fnsymbol}[1]{%
  \ensuremath{\ifcase#1\or *\or \dagger\or \ddagger\or
   \mathsection\or \|\or \bot\or \#\or @\else
   \ClassError{achemso}{Too many affiliations}
     {There are no symbols left: complain to the package
      author}\fi}}
\newcommand*{\acs@addresslist}{%
  \ifnum\acs@affilcnt>\@ne\relax
    \acs@tempcnta\@ne\relax
    \acs@address@list
  \else
    \ClassError{achemso}{No affiliations}
      {At least one affiliation is needed}%
  \fi}
\newcommand*{\acs@address@list}{%
  \advance\acs@tempcnta\@ne\relax
  \ifnum\acs@tempcnta>\acs@affilcnt\relax\else
    \acs@ifdupaffil
      {}
      {\ifnum\acs@tempcnta=\acs@affilcnt\relax
         \ifnum\acs@affilcnt>\tw@\relax
           and
         \fi
       \fi
       \csname @address@\@roman\acs@tempcnta\endcsname
       \ifnum\acs@tempcnta=\acs@affilcnt\relax\else
         ,
       \fi}%
    \acs@address@list
  \fi}
\newcommand*{\acs@ifdupaffil}{%
  \begingroup
    \acs@tempcntb\@ne\relax
    \acs@switchfalse
    \edef\acs@tempa{%
      \csname @address@\@roman\acs@tempcnta\endcsname}%
    \acs@ifdup@affil
    \expandafter\expandafter\expandafter\endgroup
    \ifacs@switch
      \expandafter\@firstoftwo
    \else
      \expandafter\@secondoftwo
    \fi}
\newcommand*{\acs@ifdup@affil}{%
  \advance\acs@tempcntb\@ne\relax
  \ifnum\acs@tempcntb=\acs@tempcnta\relax\else
    \edef\acs@tempb{%
      \csname @address@\@roman\acs@tempcntb\endcsname}%
    \ifx\acs@tempa\acs@tempb
      \expandafter\acs@switchtrue
    \fi
    \ifacs@switch\else
      \expandafter\acs@ifdup@affil
    \fi
  \fi}
\newcommand*{\acs@affilfoot}{%
  \acs@tempcnta\@ne\relax
  \acs@tempcntb\@ne\relax
  \acs@affil@foot
  \acs@tempcnta\z@\relax
  \acs@altaffil@foot}
\newcommand*{\acs@affil@foot}{%
  \advance\acs@tempcnta\@ne\relax
  \ifnum\acs@tempcnta>\acs@affilcnt\relax\else
    \advance\acs@tempcntb\@ne\relax
    \footnotetext[\acs@tempcntb]
      {\csname @affil@\@roman\acs@tempcnta\endcsname}%
    \acs@affil@foot
  \fi}
\newcommand*{\acs@altaffil@foot}{%
  \advance\acs@tempcnta\@ne\relax
  \ifnum\acs@tempcnta>\acs@authorcnt\relax\else
    \begingroup
      \@ifundefined{@altaffil@\@roman\acs@tempcnta}
        {\aftergroup\@gobble}
        {\aftergroup\@firstofone}%
    \endgroup
      {\advance\acs@tempcntb\@ne\relax
       \footnotetext[\acs@tempcntb]
         {\csname @altaffil@\@roman\acs@tempcnta\endcsname}}%
    \acs@altaffil@foot
  \fi}
\newcommand*{\acs@emaillist}{%
  \begingroup
    \renewcommand*{\UrlFont}{\sf}%
    \acs@tempcnta\z@\relax
    \acs@tempcntb\z@\relax
    \acs@email@list
    \expandafter\endgroup\expandafter\acs@tempcntb\number
      \acs@tempcntb\relax
  \ifnum\acs@tempcntb=\z@\relax
    \ClassError{achemso}{No e-mail given}
      {At lest one author must have a contact e-mail}%
  \fi}
\newcommand*{\acs@email@list}{%
  \advance\acs@tempcnta\@ne\relax
  \ifnum\acs@tempcnta>\acs@authorcnt\relax\else
    \begingroup
      \@ifundefined{@email@\@roman\acs@tempcnta}
        {\aftergroup\@gobble}
        {\aftergroup\@firstofone}%
    \endgroup
      {\advance\acs@tempcntb\@ne\relax
       \ifnum\acs@tempcntb>\@ne\relax
         ;
       \fi
       \expandafter\expandafter\expandafter\url\expandafter
         \expandafter\expandafter{%
           \csname @email@\@roman\acs@tempcnta\endcsname}}%
    \acs@email@list
  \fi}
\g@addto@macro{\document}{\maketitle}
\newfloat{scheme}{htbp}{los}
\floatname{scheme}{Scheme}
\newfloat{chart}{htbp}{loc}
\floatname{chart}{Chart}
\newfloat{graph}{htbp}{loh}
\floatname{chart}{Graph}
\newcommand*{\schemename}{Scheme}
\newcommand*{\chartname}{Chart}
\newcommand*{\graphname}{Graph}
\floatplacement{table}{htbp}
\floatplacement{figure}{htbp}
\floatstyle{plaintop}
\restylefloat{table}
\let\acs@floatboxreset\@floatboxreset
\renewcommand*{\@floatboxreset}{\centering\acs@floatboxreset}
\labelformat{scheme}{\schemename~#1}
\labelformat{chart}{\chartname~#1}
\labelformat{graph}{\graphname~#1}
\labelformat{figure}{\figurename~#1}
\labelformat{table}{\tablename~#1}
\newcommand*{\ps@achemso}{%
  \renewcommand*{\@oddfoot}{\reset@font\hfil\thepage\hfil}%
  \let\@evenfoot\@oddfoot
  \renewcommand*{\@oddhead}{%
    \reset@font
    \@author@i
    \ifnum\acs@authorcnt>\@ne\relax
      \space et al.%
    \fi
    \hfil\relax
    \ifx\@empty\@shorttitle\@empty
      \setbox\z@\hbox{\acs@title}%
      \ifdim\wd\z@>0.45\textwidth\relax
        \truncate{0.45\textwidth}{\acs@title}%
      \else
        \acs@title
      \fi
    \else
      \@shorttitle
    \fi}%
  \let\@evenhead\@oddhead}
\pagestyle{achemso}
\let\acs@startsection\@startsection
\renewcommand*{\@startsection}[6]{%
  \if@noskipsec \leavevmode \fi
  \par
  \@tempskipa #4\relax
  \@afterindenttrue
  \ifdim\@tempskipa<\z@\relax
    \@tempskipa -\@tempskipa \@afterindentfalse
  \fi
  \if@nobreak
    \everypar{}%
  \else
    \addpenalty\@secpenalty\addvspace\@tempskipa
  \fi
  \@ifstar
    {\@ssect{#3}{#4}{#5}{#6}}
    {\@ssect{#3}{#4}{#5}{#6}}}
\newcommand*{\acs@restsecnums}{%
  \let\@startsection\acs@startsection}
\let\acs@subsection\subsection
\let\acs@section\section
\newcommand*{\acs@killsecs}{%
  \newcommand*{\acs@gobblesection}{%
    \ClassWarning{achemso}
      {Sections not allowed for this manuscript type}%
    \@ifstar{\@gobble}{\@gobble}}
  \let\section\acs@gobblesection
  \let\subsection\acs@gobblesection
  \let\subsubsection\acs@gobblesection
  \AtBeginDocument{
    \renewcommand*{\bibsection}{\acs@section*{\refname}}}}
\newcommand*{\acknowledgement}{%
  \acs@subsection*{Acknowledgement}}
\newcommand*{\suppinfo}{%
  \acs@subsection*{Supporting Information Available}}
\renewcommand*{\refname}{%
  \ifnum\the\value{bibnote}>\z@\relax
    Notes and
  \fi References}
\RequirePackage{url}
\@ifundefined{acs@manuscript}
  {\newcommand*{\acs@manuscript}{article}}{}
\edef\acs@tempa{note}
\ifx\acs@manuscript\acs@tempa
  \acs@killsecs
\fi
\edef\acs@tempa{review}
\ifx\acs@manuscript\acs@tempa
  \acs@restsecnums
\fi
\edef\acs@tempa{suppinfo}
\ifx\acs@manuscript\acs@tempa
  \acs@restsecnums
  \acs@setkeys{maxauthors=0}
\fi
\if@filesw
  \acs@writebib
\fi
\endinput
%%
%% End of file `achemso.cls'.
