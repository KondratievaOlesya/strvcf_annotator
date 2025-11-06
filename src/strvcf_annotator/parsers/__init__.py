"""VCF parser module for extracting genotype and variant information."""

from .base import BaseVCFParser
from .generic import GenericParser

__all__ = ['BaseVCFParser', 'GenericParser']
