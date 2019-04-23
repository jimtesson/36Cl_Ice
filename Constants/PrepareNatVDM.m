clear variables
clc

load Nat

M0=7.74574;

Nat=Nat';
Nat(1,:)=Nat(1,:)./1000;
Nat(2,:)=Nat(2,:).*M0;